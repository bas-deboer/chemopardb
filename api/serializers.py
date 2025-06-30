from rest_framework import serializers
from protein.models import Protein, ProteinFamily
from structure.models import Structure, Entity, EntityType, StructureType, PdbData, ChemokineBindingPartner
from partner.models import Partner
from interaction.models import ChemokinePartnerIFP

# Proteins (chemokines)
class ProteinFamilySerializer(serializers.ModelSerializer):
    class Meta:
        model = ProteinFamily
        fields = ['name']
        

class ProteinSerializer(serializers.ModelSerializer):
    chemokine_ID = serializers.SerializerMethodField()
    uniprot = serializers.CharField(source='uniprot_id')
    class Meta:
        model = Protein
        fields = [
            'chemokine_ID', 'name', 'gene_name', 'subfamily', 'type',
            'species', 'full_name', 'uniprot', 'iuphar', 'sequence'
        ]

    def get_chemokine_ID(self, obj):
        return obj.id













# Structures
class StructureSerializer(serializers.ModelSerializer):
    structure_type = serializers.SlugRelatedField(
        slug_field='slug',
        queryset=StructureType.objects.all()
    )
    gene_name = serializers.SerializerMethodField()
    species = serializers.SerializerMethodField()
    pdb_code = serializers.SerializerMethodField()
    publication = serializers.SerializerMethodField()

    class Meta:
        model = Structure
        fields = ['id', 'structure_type', 'gene_name', 'species', 'resolution',
                  'publication_date', 'pdb_code', 'publication']

    def get_gene_name(self, obj):
        return obj.protein.gene_name if obj.protein else None

    def get_species(self, obj):
        return obj.protein.species if obj.protein else None

    def get_pdb_code(self, obj):
        return obj.pdb_code.index if obj.pdb_code else None

    def get_publication(self, obj):
        return obj.publication.web_link.index if obj.publication else None



# Entities
from rest_framework import serializers
from structure.models import Entity, EntityType

class EntityTypeSerializer(serializers.ModelSerializer):
    class Meta:
        model = EntityType
        fields = ['slug', 'name']

class EntitySerializer(serializers.ModelSerializer):
    entity_type = EntityTypeSerializer()  # Nested serializer to include entity_type details

    class Meta:
        model = Entity
        fields = '__all__'


class ChemokineBindingPartnerSerializer(serializers.ModelSerializer):
    class Meta:
        model = ChemokineBindingPartner
        fields = [
            'id',
            'chemokine_name',
            'chemokine_chain',
            'partner_name',
            'partner_type',
            'partner_chain',
            'partner_residues',
            'structure',
            'pdb_data',
        ]






# Partners
class PartnerSerializer(serializers.ModelSerializer):
    class Meta:
        model = Partner
        fields = '__all__'


class EntitySerializer(serializers.ModelSerializer):
    class Meta:
        model = Entity
        fields = ['id', 'name']




# Interactions
class InteractionSerializer(serializers.Serializer):
    position = serializers.IntegerField()
    type = serializers.CharField()


class ChemokinePartnerIFPSerializer(serializers.ModelSerializer):
    structure_id = serializers.IntegerField(source='structure.id', read_only=True)
    binding_pair_id = serializers.IntegerField(source='binding_pair.id', read_only=True)

    class Meta:
        model = ChemokinePartnerIFP
        fields = ['structure_id', 'binding_pair_id', 'ifp_string']






class StructureTypeSerializer(serializers.ModelSerializer):
    class Meta:
        model = StructureType
        fields = ['slug']


class PdbDataSerializer(serializers.ModelSerializer):
    class Meta:
        model = PdbData
        fields = '__all__'



class ErrorResponseSerializer(serializers.Serializer):
    code = serializers.IntegerField()
    message = serializers.CharField()