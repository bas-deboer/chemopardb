import requests
import logging
from django.db import models
from django.db.utils import IntegrityError
from django.utils.text import slugify

logger = logging.getLogger(__name__)

class Publication(models.Model):
    web_link = models.OneToOneField('common.WebLink', null=True, on_delete=models.CASCADE)
    journal = models.ForeignKey('PublicationJournal', on_delete=models.CASCADE, null=True)
    title = models.TextField(null=True)
    authors = models.TextField(null=True)
    year = models.IntegerField(null=True)
    reference = models.TextField(null=True)

    def __str__(self):
        return "{!s} ({!s})".format(self.journal, self.year)

    class Meta:
        db_table = 'publication'

    @staticmethod
    def get_or_create_from_type(identifier, wr):
        if isinstance(identifier, int):
            identifier = str(identifier)
        if len(identifier) > 0:
            try:
                wl = WebLink.objects.get(index__iexact=identifier, web_resource=wr)
                publications = Publication.objects.filter(web_link=wl)
                if publications.count() == 1:
                    return publications.first()
            except WebLink.DoesNotExist:
                wl = False

            pub = Publication()
            pub.web_link, created = WebLink.objects.get_or_create(defaults={"index": identifier}, index__iexact=identifier, web_resource=wr)
            if wr.slug == "doi":
                pub.update_from_doi(identifier)
            elif wr.slug == "pubmed":
                pub.update_from_pubmed_data(identifier)
            pub.save()
            return pub
        else:
            return False

    @classmethod
    def get_or_create_from_doi(cls, doi):
        return cls.get_or_create_from_type(doi, WebResource.objects.get(slug="doi"))

    @classmethod
    def get_or_create_from_pubmed(cls, pmid):
        return cls.get_or_create_from_type(pmid, WebResource.objects.get(slug="pubmed"))

    @staticmethod
    def update_from_doi(doi):
        url = f'http://api.crossref.org/works/{doi}'
        try:
            response = requests.get(url)
            response.raise_for_status()
            pub = response.json()
        except requests.exceptions.RequestException as e:
            logger.error(f"Error fetching data from CrossRef: {e}")
            return False

        if pub:
            try:
                title = pub['message']['title'][0]
                try:
                    year = pub['message']['created']['date-parts'][0][0]
                except KeyError:
                    year = pub['message'].get('deposited', {}).get('date-parts', [[None]])[0][0]

                authors = ['{} {}'.format(x['family'], ''.join([y[:1] for y in x['given'].split()]))
                           for x in pub['message']['author']]
                authors = ', '.join(authors)

                reference = {}
                for f in ['volume', 'page']:
                    reference[f] = pub['message'].get(f, 'X')
                reference = '{}:{}'.format(reference['volume'], reference['page'])

                journal = pub['message']['container-title'][0]
                journal_abbr = pub['message'].get('short-container-title', [slugify(journal)])[0]

                try:
                    journal, created = PublicationJournal.objects.get_or_create(defaults={"name": journal, 'slug': journal_abbr}, name__iexact=journal)
                    if created:
                        logger.info('Created journal {}'.format(journal))
                except IntegrityError:
                    journal = PublicationJournal.objects.get(name=journal)
            except Exception as msg:
                logger.warning('Processing data from CrossRef for {} failed: {}'.format(doi, msg))
        else:
            print("Publication not on crossref", doi)
        return False

class PublicationJournal(models.Model):
    slug = models.CharField(null=True, blank=True, max_length=255)
    name = models.TextField(unique=True)

    def __str__(self):
        return "{!s} ({!s})".format(self.name, self.slug)

    class Meta:
        db_table = 'publication_journal'

class WebLink(models.Model):
    web_resource = models.ForeignKey('WebResource', on_delete=models.CASCADE)
    index = models.TextField()

    def __str__(self):
        return f"{self.index}"

    class Meta:
        db_table = 'web_link'

class WebResource(models.Model):
    slug = models.SlugField()
    name = models.TextField()
    url = models.TextField()

    def __str__(self):
        return '<{}>'.format(self.name)

    class Meta:
        db_table = 'web_resource'

class ResiduePosition(models.Model):
    position = models.PositiveIntegerField(unique=True)
    ccn_number = models.CharField(max_length=32, unique=True)
    
    def __str__(self):
        return f"{self.number}: {self.ccn_number}"