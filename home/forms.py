from django import forms
from django.core.validators import EmailValidator

class DateForm(forms.Form):
    start = forms.DateField(widget=forms.DateInput(attrs={'type': 'date'}))
    end = forms.DateField(widget=forms.DateInput(attrs={'type': 'date'}))


class ChemokineDiagramForm(forms.Form):
    """Form allowing users to paste a sequence that will be segmented interactively."""

    sequence = forms.CharField(
        widget=forms.Textarea(attrs={"rows": 5, "placeholder": "Enter protein sequence"}),
        label="Sequence",
    )
