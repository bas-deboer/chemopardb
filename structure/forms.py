from django import forms

class ChainSelectionForm(forms.Form):
    chain = forms.ChoiceField(choices=[('A', 'Chain A'), ('B', 'Chain B'), ('C', 'Chain C'), ('D', 'Chain D')], label="Select Chain")
