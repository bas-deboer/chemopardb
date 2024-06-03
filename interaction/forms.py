from django import forms

class PDBSelectionForm(forms.Form):
    pdb_ids = forms.CharField(label='Enter PDB IDs', widget=forms.TextInput(attrs={'placeholder': 'e.g., 4HHB, 1GZX'}))
