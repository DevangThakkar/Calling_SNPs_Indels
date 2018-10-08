from django import forms

import sys
sys.path.append('/home/vcm/')

import script_cbb520

class ContactForm(forms.Form):
    OPTIONS_DATA = (
            ("a", "Use existing data"),
            ("b", "Download and process the data all over again (Takes time)     "),
            )
    data = forms.ChoiceField(widget=forms.Select(choices=OPTIONS_DATA),
                                        choices=OPTIONS_DATA)

    OPTIONS_PRINT = (
            ("a", "Print number of high quality pairs (Takes ~10 minutes)     "),
            ("b", "Print number of SNPs"),
            ("c", "Print types of SNPs"),
            ("d", "Print number of single-base Indels"),
            ("e", "Print types of Indels"),
            ("f", "Print mean/std dev of SNPs per 10kb window"),
            ("g", "Print mean/std dev of Indels per 10kb window"),
            ("h", "Print expected vs observed number of windows (SNPs)"),
            ("i", "Print expected vs observed number of windows (Indels)")
            )
    options = forms.MultipleChoiceField(widget=forms.CheckboxSelectMultiple,
                                        choices=OPTIONS_PRINT)

    def clean(self):
        cleaned_data = super(ContactForm, self).clean()
        data = cleaned_data.get('data')
        options = cleaned_data.get('options')
        # call scripts from here
        name = 'SRR4841864'
        index = 'sacCer3bwaidx'
        
        return_str = ''
        if 'b' in data:
            return_str += script_cbb520.delete_files(name, index)
            return_str += script_cbb520.process(name, index)

        if 'a' in options:
            return_str += script_cbb520.count_high_qual(name, 50, 25)
        if 'b' in options:
            return_str += script_cbb520.calc_snp_types(name, print_num=True)
        if 'c' in options:
            return_str += script_cbb520.calc_snp_types(name, print_type=True)
        if 'd' in options:
            return_str += script_cbb520.calc_indel_types(name, print_num=True)
        if 'e' in options:
            return_str += script_cbb520.calc_indel_types(name, print_type=True)
        if 'f' in options:
            return_str += script_cbb520.calc_snp_types(name, print_stats=True)
        if 'g' in options:
            return_str += script_cbb520.calc_indel_types(name, print_stats=True)
        if 'h' in options:
            return_str += script_cbb520.calc_snp_types(name, print_exp_obs=True)
        if 'i' in options:
            return_str += script_cbb520.calc_indel_types(name, print_exp_obs=True)
        return return_str