# -*- coding: utf-8 -*-
"""
Created on Mon Oct  8 14:16:25 2012

@author: proto

"""

import networkx as nx
import json

import bioservices
import pyparsing as pyp

from twisted.web import xmlrpc, server
from twisted.internet import threads,reactor


import concurrent.futures
import urllib2

import tempfile
import subprocess
from os import getcwd,remove
# Restrict to a particular path.
import datetime

# Create server
port = 9200
bngDistro  = '/home/ubuntu/bionetgen/bng2/'
#bngDistro = '/home/proto/workspace/bionetgen/bng2/'

#server = SimpleXMLRPCServer(("127.0.0.1", port), requestHandler=RequestHandler)


#def updateDictionary(tmpArray,key,newInfo):

def resolveAnnotations(annotations):
    tmpArray = {}
    futures = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=8) as executor:
        for annotation in annotations:
            futures.append(executor.submit(resolveAnnotation, annotation))
        for future in concurrent.futures.as_completed(futures):
            resolvedAnnotation = future.result()
    
            for element in resolvedAnnotation:
                if resolvedAnnotation[element] != annotation:
                    if element not in tmpArray:
                        tmpArray[element] = []
                    if isinstance(resolvedAnnotation[element],list):
                        tmpArray[element].extend(resolvedAnnotation[element])
                    else:
                        tmpArray[element].append(resolvedAnnotation[element])
    finalArray = []
    #transforming to an [key [result],....] structure since the ouput cant be dictionaries
    for element in tmpArray:
        finalArray.append([element,tmpArray[element]])
    print finalArray
    return finalArray

def removeTags(taggedInformation):
        taggedInformation2 = taggedInformation.decode('ascii','ignore')
        try:
            goGrammar = pyp.Suppress('<' + pyp.Word(pyp.alphanums) + '>') +  pyp.Word(pyp.alphanums + pyp.alphas8bit + ' .,_-') + pyp.Suppress('</' + pyp.Word(pyp.alphanums) + '>')
            tmp = goGrammar.parseString(str(taggedInformation2))
        except pyp.ParseException:
            tmp = [taggedInformation2]
        return tmp[0]
    
def resolveAnnotation(annotation):
    if not hasattr(resolveAnnotation, 'db'):
        servicesSignature = {'ch':bioservices.ChEBI,'bio':bioservices.BioModels,'uni':bioservices.UniProt,
                                                        'k':bioservices.kegg.KeggParser,'qg':bioservices.QuickGO}

        resolveAnnotation.services = {}
        for signature in servicesSignature:
            try:
                resolveAnnotation.services[signature] = servicesSignature[signature](verbose=False)
            except:
                continue
        resolveAnnotation.db = {}
        '''
        resolveAnnotation.ch = bioservices.ChEBI(verbose=False)
        resolveAnnotation.bio = bioservices.BioModels(verbose=False)
        resolveAnnotation.uni = bioservices.UniProt(verbose=False)
        resolveAnnotation.k = bioservices.kegg.KeggParser(verbose=False)
        resolveAnnotation.qg = bioservices.QuickGO(verbose=False)
        '''
        resolveAnnotation.s = bioservices.Service('name')
 

        resolveAnnotation.db['http://identifiers.org/uniprot/P62988'] = 'http://identifiers.org/uniprot/P62988'
        resolveAnnotation.db['http://identifiers.org/uniprot/P06842'] = 'http://identifiers.org/uniprot/P06842'
        resolveAnnotation.db['http://identifiers.org/uniprot/P07006'] = 'http://identifiers.org/uniprot/P06842'
        
    if annotation in resolveAnnotation.db:
        return resolveAnnotation.db[annotation]
    
        
    tAnnotation = annotation.replace('%3A',':')
    tAnnotation = annotation.split('/')[-1]
    #tAnnotation = re.search(':([^:]+:[^:]+$)',tAnnotation).group(1)
    tmpArray = {}
    try:
        if 'obo.go' in annotation and 'qg' in resolveAnnotation.services:
            res = resolveAnnotation.services['qg'].Term(tAnnotation)
            tmp = res.findAll('name')
            finalArray = []
            goGrammar = pyp.Suppress(pyp.Literal('<name>')) +  pyp.Word(pyp.alphanums + ' -_') + pyp.Suppress(pyp.Literal('</name>')) 
            for x in tmp:
                try:
                    tagString = str(goGrammar.parseString(str(x))[0])
                    if tagString not in ['Systematic synonym']:
                        finalArray.append(str(goGrammar.parseString(str(x))[0]))
                except pyp.ParseBaseException:
                    continue
            tmp = finalArray
            resolveAnnotation.db[annotation] = tmp
            tmpArray['tags'] =  tmp

        elif 'kegg' in annotation and 'k' in resolveAnnotation.services:
            
            data = resolveAnnotation.k.get(tAnnotation)
            dict_data =  resolveAnnotation.services['k'].parse(data)
            resolveAnnotation.db[annotation] = dict_data['name']
            tmpArray['tags'] =  resolveAnnotation.db[annotation]
            
        elif 'uniprot' in annotation and 'uni' in resolveAnnotation.services:
            identifier = annotation.split('/')[-1]
            result = resolveAnnotation.services['uni'].quick_search(identifier)
            resolveAnnotation.db[annotation] = result[result.keys()[0]]['Protein names']
            tmpArray['tags'] = resolveAnnotation.db[annotation]
            
        elif 'chebi' in annotation and 'ch' in resolveAnnotation.services:
            tmp = annotation.split('/')[-1]
            
            entry = resolveAnnotation.services['ch'].getLiteEntity(tmp)
            for element in entry:
                resolveAnnotation.db[annotation] = element['chebiAsciiName']
                tmpArray['tags'] =  resolveAnnotation.db[annotation]
        elif 'cco' in annotation or 'pirsf' in annotation or 'pubchem' in annotation or 'omim' in annotation:
            tmpArray['ignore'] =  annotation
        elif 'biomodels.db' in annotation and 'bio' in resolveAnnotation.services:
            tmp = annotation.split('/')[-1]
            tmpArray = {}
            try:
                request = resolveAnnotation.services['bio'].getSimpleModelsByIds(tmp)
            except ValueError:
                tmpArray['ignore'] = annotation
                return tmpArray
            entry = resolveAnnotation.s.easyXML(request)
            tmpArray['modelName'] = []
            for tag in entry['modelname']:
                tmpArray['modelName'].append(removeTags(str(tag)))
            tmpArray['author'] = []
            for tag in entry['author']:
                tmpArray['author'].append(removeTags(str(tag)))
            
        #elif 'taxonomy' in annotation:
            #uniprot stuff for taxonomy
        #    pass
            '''
            url = 'http://www.uniprot.org/taxonomy/'
            params = {
            'from':'ACC',
            'to':'P_REFSEQ_AC',
            'format':'tab',
            'query':'P13368 P20806 Q9UM73 P97793 Q17192'
            }
            
            data = urllib.urlencode(params)
            request = urllib2.Request(url, data)
            contact = "" # Please set your email address here to help us debug in case of problems.
            request.add_header('User-Agent', 'Python contact')
            response = urllib2.urlopen(request)
            page = response.read(200000)
            '''
        else:
            print 'ERRROERROROERRRO',annotation
            #assert(False)
            tmpArray['ignore'] =  annotation
    except (SyntaxError,TypeError,urllib2.HTTPError):
        tmpArray['ignore'] = annotation
        #finally:
    return tmpArray

def generateTimeSeries(bnglFile):
    import signal
    import os
    import time
    pointer = tempfile.mkstemp(suffix='.bngl',text=True)

    timeout = 500
    name = pointer[1].split('.')[0]

    with open(pointer[1],'w' ) as f:
        f.write(bnglFile)
    try:
        
        start = datetime.datetime.now()
        print start
        result = subprocess.Popen(['perl',bngDistro +'BNG2.pl','--outdir={0}'.format(os.sep.join(pointer[1].split(os.sep)[0:-1]) + os.sep),pointer[1]])
        while result.poll() is None:
            time.sleep(0.1)
            now = datetime.datetime.now()
            if (now - start).seconds > timeout:
                os.kill(result.pid, signal.SIGKILL)
                os.waitpid(-1, os.WNOHANG)
                raise OSError
    except OSError:
            return None
    fileName = '{0}.gdat'.format(name)
    return fileName

def generateContactMap(bnglFile,graphType):
    import signal
    import os
    import time
    pointer = tempfile.mkstemp(suffix='.bngl',text=True)
    name = pointer[1].split('.')[0]
    timeout = 120
    with open(pointer[1],'w' ) as f:
        f.write(bnglFile)
    try:
        
        start = datetime.datetime.now()
        print start
        result = subprocess.Popen(['perl',bngDistro +'Perl2/visualize.pl',pointer[1],
                graphType,graphType],cwd=bngDistro+'Perl2/')
        while result.poll() is None:
            time.sleep(0.1)
            now = datetime.datetime.now()
            if (now - start).seconds > timeout:
                os.kill(result.pid, signal.SIGKILL)
                os.waitpid(-1, os.WNOHANG)
                raise OSError
        
        #subprocess.call(['perl',bngDistro +'visualize.pl',pointer[1],
        #       graphType,graphType],cwd=bngDistro)
    except OSError:
        #print now
        #TODO: we have to return a proper error message. Right now it's just an empty file
        #alternatively we could recognize empty files as error messages
        print 'system did not finish. Returning empty file'
        with open('{1}_{0}.gml'.format(graphType,name),'w') as f:
            f.write('graph\n[]')

    #finally:
    #    remove(pointer[1])
    fileName =    '{1}_{0}.gml'.format(graphType,name)
    return fileName


def inheritColor(gmlText,node):
    pass


def gml2cyjson(gmlText):
    import random
    r = lambda: random.randint(0,255)

    jsonDict = {}
    jsonDict['style'] =  [
    {
      'selector': 'node',
      'css': {
        'content': 'data(label)',
        'text-valign': 'center',
        'text-halign': 'center'
      }
    },
    {
      'selector': '$node > node',
      'css': {
        'padding-top': '20px',
        'padding-left': '10px',
        'padding-bottom': '10px',
        'padding-right': '10px',
        'text-valign': 'top',
        'text-halign': 'center'
      }
    },
    {
      'selector': 'edge',
      'css': {
        'target-arrow-shape': 'triangle'
      }
    },
    {
      'selector': ':selected',
      'css': {
        'background-color': 'black',
        'line-color': 'black',
        'target-arrow-color': 'black',
        'source-arrow-color': 'black'
      }
    }
  ]
    jsonDict['elements'] = {'nodes':[]}
    jsonDict['elements']['edges'] = []
    #for nd in gmlText.node:
    #    gmlText.node[nd].pop('id')
    #jsgrph = json_graph.node_link_data(gmlText)
    colorDict = {}
    for node in gmlText.node:
        if gmlText.node[node] == {}:
            continue
        tmp = {'data':{}}
        tmp['data']['id'] = str(node)
        tmp['data']['label'] = str(gmlText.node[node]['label'])
        
        if 'gid' in gmlText.node[node]:
            tmp['data']['parent'] =  str(gmlText.node[node]['gid'])
            if str(gmlText.node[node]['gid']) not in colorDict:
                if 'gid' in gmlText.node[str(gmlText.node[node]['gid'])]:
                    if str(gmlText.node[str(gmlText.node[node]['gid'])]['gid']) not in colorDict:
                        newColor = '#%02X%02X%02X' % (r(),r(),r())
                        colorDict[str(gmlText.node[str(gmlText.node[node]['gid'])]['gid'])] = newColor
                        colorDict[str(gmlText.node[node]['gid'])] = newColor
                    else:
                        colorDict[str(gmlText.node[node]['gid'])] = colorDict[str(gmlText.node[str(gmlText.node[node]['gid'])]['gid'])] 
                else:
                    colorDict[str(gmlText.node[node]['gid'])] = '#%02X%02X%02X' % (r(),r(),r())
            colorDict[str(node)] = colorDict[str(gmlText.node[node]['gid'])]
        if str(node) not in colorDict:
            colorDict[str(node)] = '#%02X%02X%02X' % (r(),r(),r())
        tmp['data']['faveColor'] = colorDict[str(node)]

        jsonDict['elements']['nodes'].append(tmp)
    for link in gmlText.edge:
        for dlink in gmlText.edge[link]:
            if link != '' and dlink != '':
                tmp = {'data':{}}
                tmp['data']['source'] = int(link)
                tmp['data']['target'] = int(dlink)
                tmp['data']['id'] = '{0}_{1}'.format(tmp['data']['source'],tmp['data']['target'])
                tmp['data']['faveColor'] = colorDict[str(link)]
                jsonDict['elements']['edges'].append(tmp)

    jsonDict['layout'] = {
    'name': 'cose',
    'padding': 4,
     'fit'                 : True, 
   'nodeRepulsion'       : 10000,
    'nodeOverlap'         : 10,
    'idealEdgeLength'     : 10,
   'edgeElasticity'      : 100,
    'nestingFactor'       : 5, 
    'gravity'             : 250, 
    'numIter'             : 100,
    'initialTemp '        : 200,
    'coolingFactor'       : 0.95, 
    'minTemp'             : 1  },
  
    jsonDict['ready'] =  'function(){window.cy = this;}'

    return jsonDict

def gdat2c3json(gdatText):
    #from collections import OrderedDict
    gdatlines = gdatText.split('\n')
    labels = gdatlines[0].split()[1:]
    dataPoints = [x.split() for x in gdatlines[1:]]
    c3json = {'xs':{},'columns':[]}
    for idx,label in enumerate(labels):
        tmp = [label]
        for dataPoint in dataPoints:
            if dataPoint == []:
                continue
            tmp.append(dataPoint[idx])
        c3json['columns'].append(tmp)
    for label in labels[1:]:
        c3json['xs'][label] = labels[0]
    return c3json
     

def tmpGenerateCont(bnglFile,graphType):
        fileName = generateContactMap(bnglFile, graphType)
        try:
            gml = nx.read_gml(fileName)
            result = gml2cyjson(gml)
        except IOError:
            result = '{}'
        
        with open('out.json','w') as f:
            return json.dumps(result,f,indent=1, separators=(',', ': '))
        

def simulateFile(bnglFile):
    pass

def getContactMap(bnglFile,graphType):
        '''
        calls John's bgn map generator
        '''
        fileName = generateContactMap(bnglFile, graphType)
        try:
            with open(fileName) as f:
                gmlText = f.read()
            gml = nx.read_gml(fileName)
        except IOError:
            return {'jsonStr':'','gmlStr':''}
        result = gml2cyjson(gml)
        jsonStr = json.dumps(result,indent=1, separators=(',', ': '))
        result = {'jsonStr':jsonStr,'gmlStr':gmlText}
        print gmlText
        #remove(fileName)
        return result
    
def getTimeSeries(bnglFile):
        fileName = generateTimeSeries(bnglFile)
        try:
            if fileName == None:
                raise IOError
            with open(fileName) as f:
                gdatText = f.read()
        except IOError:
            return {'jsonStr':'','gdatStr':''}
        result = gdat2c3json(gdatText)
        jsonStr = json.dumps(result,indent=1,separators=(',',':'))
        result = {'jsonStr':jsonStr,'gdatStr':gdatText}
        return result
class AnnotationServer(xmlrpc.XMLRPC):


    def xmlrpc_resolveAnnotations(self,annotations):

        #result = threads.deferToThread(resolveAnnotations,annotations)
        result = resolveAnnotations(annotations)
        return result
  
    def xmlrpc_getContactMap(self,bnglFile,graphType):
        #result = threads.deferToThread(getContactMap,bnglFile,graphType)
        result = getContactMap(bnglFile,graphType)
        return result
    
    def xmlrpc_getTimeSeries(self,bnglFile):
        result = getTimeSeries(bnglFile)
        return result

#server.register_function(is_even, "is_even")

#resolveAnnotations(tmpD)        
if __name__ == '__main__':
    print "Listening on port {0}...".format(port)
    r = AnnotationServer()
    reactor.listenTCP(port, server.Site(r))
    reactor.run()
    #print getcwd()
    #gml = nx.read_gml('/tmp/tmpy0ug0r_contact.gml')

    #gml2cyjson(gml) 
    #a = ['http://identifiers.org/pubmed/10514507', 'http://identifiers.org/biomodels.db/BIOMD0000000048', 'http://identifiers.org/biomodels.db/MODEL6624193277', 'http://identifiers.org/taxonomy/10116']
    #for element in a:
    #   print element
    #   print resolveAnnotation(element)

        

