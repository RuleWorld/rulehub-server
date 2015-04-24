from google.appengine.ext import ndb

class PublicationInfo(ndb.Model):
    name = ndb.StringProperty()
    journal = ndb.StringProperty()

class AnnotationInfo(ndb.Model):
    database = ndb.StringProperty()
    databaseID = ndb.StringProperty()

class ContentInfo(ndb.Model):
    content = ndb.BlobKeyProperty() #BlobInfo(blobkey)
    date = ndb.DateTimeProperty(auto_now_add=True)
    contactMap = ndb.BlobKeyProperty()
    contactMapJson = ndb.JsonProperty()
    processMap = ndb.BlobKeyProperty()
    timeSeries = ndb.BlobKeyProperty()
    timeSeriesJson = ndb.JsonProperty()
    processMapJson = ndb.JsonProperty()     

    @classmethod
    def create(cls, params, doc_id):
        content = cls(content=params['content'])
        if 'timeSeries' in params:
            content.timeSeries = params['timeSeries']
            content.timeSeriesJson = params['timeSeriesJson']
        if 'processMap' in params:
            content.processMap=params['processMap']
            content.processMapJson=params['processMapJson']
        if 'contactMap' in params:
            content.contactMap=params['contactMap']
            content.contactMapJson=params['contactMapJson']
        return content


class ModelInfo(ndb.Model):
    """Models an individual Guestbook entry with author, content, and date."""
    author = ndb.StringProperty(repeated=True)
    '''
    content = ndb.BlobKeyProperty() #BlobInfo(blobkey)
    contactMap = ndb.BlobKeyProperty()
    contactMapJson = ndb.JsonProperty()
    processMap = ndb.BlobKeyProperty()
    timeSeries = ndb.BlobKeyProperty()
    timeSeriesJson = ndb.JsonProperty()
    processMapJson = ndb.JsonProperty()     
    '''
    contentHistory = ndb.StructuredProperty(ContentInfo,repeated=True)
    name = ndb.StringProperty()
    description = ndb.StringProperty()
    date = ndb.DateTimeProperty(auto_now_add=True)
    publication = ndb.StructuredProperty(PublicationInfo)
    fileFormat = ndb.StringProperty(choices=set(["bngl","kappa"]))
    submitter = ndb.StringProperty()
    annotationInfo = ndb.KeyProperty(kind=AnnotationInfo,repeated=True)
    tags = ndb.StringProperty(repeated=True)
    structuredTags = ndb.StringProperty(repeated=True)
    privacy = ndb.StringProperty()
    notes = ndb.StringProperty()

    doc_id = ndb.StringProperty()

    def update_core(self, params, doc_id):
        """Update 'core' values from the given params dict and doc_id."""
        self.populate(doc_id=doc_id)

    @classmethod
    def create(cls, params, doc_id):
        """Create a new product entity from a subset of the given params dict
        values, and the given doc_id."""
        prod = cls(id=doc_id,
            content=params['content'],  doc_id=doc_id,name=params['name'],
            submitter=params['submitter'],privacy=params['privacy'],tags=params['tags'],
            structuredTags=params['structuredTags'],author=params['author']
            )
        if 'notes' in params:
            prod.notes = params['notes']
        if 'tags' in params:
            prod.tags=params['tags']
        content = ContentInfo.create(params)
        prod.contentHistory = [content]
        #if 'structuredTags' in params:
         #   prod.structuredTags=params['structuredTags'],
        #prod.put()
        return prod

