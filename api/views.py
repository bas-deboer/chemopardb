from django.shortcuts import render
from django.http import JsonResponse
from django.db.models import Q
from drf_spectacular.utils import extend_schema, OpenApiResponse, OpenApiParameter
from rest_framework import viewsets
from rest_framework.response import Response
from rest_framework.renderers import JSONRenderer
from rest_framework import generics
from rest_framework.views import APIView

from structure.models import Structure, Entity
from protein.models import Protein, ProteinFamily
from partner.models import Partner
from api.serializers import StructureSerializer, EntitySerializer, ProteinFamilySerializer, ProteinSerializer, PartnerSerializer, ErrorResponseSerializer, InteractionSerializer
from interaction.models import ChemokinePartnerIFP
from api.serializers import ChemokinePartnerIFPSerializer
from api.serializers import EntitySerializer
from structure.models import ChemokineBindingPartner
from api.serializers import ChemokineBindingPartnerSerializer

# Autocomplete search bar
def protein_autocomplete(request):
    query = request.GET.get('term', '')
    proteins = Protein.objects.filter(gene_name__icontains=query)[:20]
    results = [{'label': protein.gene_name, 'type': 'protein', 'url': f'/protein/{protein.gene_name}/'} for protein in proteins]
    return JsonResponse(results, safe=False)

def structure_autocomplete(request):
    query = request.GET.get('term', '')
    structures = Structure.objects.filter(
        Q(pdb_code__index__icontains=query) | Q(protein__gene_name__icontains=query))[:20]
    results = [{'label': f"{structure.pdb_code.index} ({structure.protein.gene_name})", 'type': 'structure', 'url': f'/structure/{structure.id}/'} for structure in structures]
    return JsonResponse(results, safe=False)
    
    
@extend_schema(
    tags=['Proteins'],
    summary="List all protein details",
    description="Returns a list of all available protein details stored in ChemoPar-db.",
    responses={
        200: OpenApiResponse(
            response={
                'type': 'array',
                'items': {'type': 'object'}
            },
            description="A list of protein details"
        )
    }
)
class ProteinDetailList(generics.ListAPIView):
    queryset = Protein.objects.all()
    serializer_class = ProteinSerializer


# Structures
@extend_schema(
    tags=['Structures'],
    summary="Get a list of all structures or details of a specific structure",
    description="The Structure list endpoint returns a list of all chemokine structures or details of a specific structure if an ID is provided.",
    parameters=[
        OpenApiParameter(name='id', description='ID of the structure to retrieve', required=False, type=int),
    ],
    responses={
        200: OpenApiResponse(
            response=StructureSerializer(many=True),
            description="A list of structures or a single structure based on the provided ID"
        ),
        'default': OpenApiResponse(
            response=ErrorResponseSerializer,
            description="Unexpected error"
        )
    }
)
class StructureDetailList(generics.ListAPIView):
    """
    Returns a list of all available structures or details of a specific structure if an ID is provided.
    """
    queryset = Structure.objects.all()
    serializer_class = StructureSerializer

    def get_queryset(self):
        queryset = super().get_queryset()
        structure_id = self.request.query_params.get('id', None)
        if structure_id:
            queryset = queryset.filter(id=structure_id)
        return queryset



@extend_schema(
    tags=['Structures'],
    summary="List all structures",
    description="Get a list of all structures.",
    responses={
        200: OpenApiResponse(
            response=StructureSerializer(many=True),
            description="A list of structures"
        )
    }
)
class StructureList(generics.ListAPIView):
    queryset = Structure.objects.all()
    serializer_class = StructureSerializer


@extend_schema(
    tags=['Structures'],
    summary="Get structure details",
    description="Retrieve details of a single structure by its ID.",
    responses={
        200: OpenApiResponse(
            response=StructureSerializer,
            description="Details of a single structure"
        )
    }
)
class StructureDetail(generics.RetrieveAPIView):
    queryset = Structure.objects.all()
    serializer_class = StructureSerializer
    lookup_field = 'pk'


@extend_schema(
    tags=['Structures'],
    summary="Search structures",
    description="Search for structures by PDB code.",
    responses={
        200: OpenApiResponse(
            response=StructureSerializer(many=True),
            description="A list of structures"
        )
    }
)
class StructureSearch(generics.ListAPIView):
    serializer_class = StructureSerializer

    def get_queryset(self):
        queryset = Structure.objects.all()
        pdb_code = self.request.query_params.get('pdb_code', None)
        if pdb_code:
            queryset = queryset.filter(pdb_code__icontains=pdb_code)
        return queryset



# Binding pairs
@extend_schema(
    tags=['Binding Pairs'],
    summary="List all binding pairs",
    description="Get a list of all binding pairs.",
    responses={
        200: OpenApiResponse(
            response=EntitySerializer(many=True),
            description="A list of binding pairs"
        )
    }
)
class EntityListAPIView(generics.ListAPIView):
    """
    API endpoint to list all binding pairs (entities) or filter by specific parameters.
    """
    queryset = Entity.objects.all()
    serializer_class = EntitySerializer

    def get_queryset(self):
        queryset = super().get_queryset()
        # Example: Filtering by name, organism, or other parameters
        name = self.request.query_params.get('name', None)
        if name:
            queryset = queryset.filter(name__icontains=name)
        return queryset



@extend_schema(
    tags=['Binding Pairs'],
    summary="List all chemokine binding partner pairs",
    description="Returns a list of all chemokine binding partner pairs.",
    responses={
        200: OpenApiResponse(
            response=ChemokineBindingPartnerSerializer(many=True),
            description="A list of chemokine binding partner pairs"
        ),
        'default': OpenApiResponse(
            response=ErrorResponseSerializer,
            description="Unexpected error"
        )
    }
)
class ChemokineBindingPartnerList(generics.ListAPIView):
    queryset = ChemokineBindingPartner.objects.all()
    serializer_class = ChemokineBindingPartnerSerializer





# Interactions
class InteractionListView(APIView):
    @extend_schema(
        tags=['Interactions'],
        summary="Get possible interactions",
        description="Returns a list of possible interactions with their positions and types.",
        responses={
            200: OpenApiResponse(
                response=InteractionSerializer(many=True),
                description="A list of possible interactions"
            )
        }
    )
    def get(self, request):
        possible_interactions = [
            "Hydrophobic", "Pistacking", "HBDonor", "HBAcceptor",
            "Anionic", "Cationic", "CationPi", "PiCation", "VdWContact"
        ]
        interaction_list = [
            {"position": index + 1, "type": interaction}
            for index, interaction in enumerate(possible_interactions)
        ]
        return Response(interaction_list)
    
    

@extend_schema(
    tags=['Interactions'],
    summary="List chemokine partner IFPs or filter by structure ID",
    description="Returns a list of all chemokine partner IFPs or filters by a specific structure ID if provided.",
    parameters=[
        OpenApiParameter(
            name='structure_id',
            description='ID of the structure to filter by',
            required=False,
            type=int
        ),
    ],
    responses={
        200: OpenApiResponse(
            response=ChemokinePartnerIFPSerializer(many=True),
            description="A list of chemokine partner IFPs or a single IFP based on the provided structure ID"
        ),
        'default': OpenApiResponse(
            response=ErrorResponseSerializer,
            description="Unexpected error"
        )
    }
)
class ChemokinePartnerIFPList(generics.ListAPIView):
    queryset = ChemokinePartnerIFP.objects.all()
    serializer_class = ChemokinePartnerIFPSerializer

    def get_queryset(self):
        queryset = super().get_queryset()
        structure_id = self.request.query_params.get('structure_id', None)
        if structure_id:
            queryset = queryset.filter(structure__id=structure_id)
        return queryset







# Entities
@extend_schema(
    tags=['Entities'],
    summary="List all entities",
    description="Get a list of all entities.",
    responses={
        200: OpenApiResponse(
            response=EntitySerializer(many=True),
            description="A list of entities"
        )
    }
)
class EntityList(generics.ListAPIView):
    queryset = Entity.objects.all()
    serializer_class = EntitySerializer


@extend_schema(
    tags=['Entities'],
    summary="Get entity details",
    description="Retrieve details of a single entity by its ID.",
    responses={
        200: OpenApiResponse(
            response=EntitySerializer,
            description="Details of a single entity"
        )
    }
)
class EntityDetail(generics.RetrieveAPIView):
    queryset = Entity.objects.all()
    serializer_class = EntitySerializer










# Partners
@extend_schema(
    tags=['Partners'],
    summary="Get all partners",
    description="Returns a list of all available partners.",
    responses={
        200: OpenApiResponse(
            response={
                'type': 'array',
                'items': {'type': 'object'}
            },
            description="A list of partner details"
        )
    }
)
class PartnersList(generics.ListAPIView):
    queryset = Partner.objects.all()
    serializer_class = PartnerSerializer


#@extend_schema(
#    tags=['Partners'],
#    summary="List all chemokine partner pairs or filter by ID",
#    description="Returns a list of all chemokine partner pairs or filters by a specific chemokine partner pair ID if provided.",
#    parameters=[
#        OpenApiParameter(name='chemokinepartnerpair_ID', description='ID of the chemokine partner pair to retrieve', required=False, type=int),
#    ],
#    responses={
#        200: OpenApiResponse(
#            response=ChemokinePartnerPairSerializer(many=True),
#            description="A list of chemokine partner pairs or a single pair based on the provided ID"
#        ),
#        'default': OpenApiResponse(
#            response=ErrorResponseSerializer,
#            description="Unexpected error"
#        )
#    }
#)
#class ChemokinePartnerPairList(generics.ListAPIView):
#    """
#    Returns a list of all chemokine partner pairs or details of a specific pair if an ID is provided.
#    """
#    queryset = ChemokinePartnerPair.objects.all()
#    serializer_class = ChemokinePartnerPairSerializer
#
#    def get_queryset(self):
#        queryset = super().get_queryset()
#        pair_id = self.request.query_params.get('chemokinepartnerpair_ID', None)
#        if pair_id:
#            queryset = queryset.filter(id=pair_id)
#        return queryset