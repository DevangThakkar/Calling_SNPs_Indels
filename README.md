# CBB520 Assignment
## Author: Devang Thakkar

The repository here contains the code written for the assignment. The python script script_cbb520.py contains all the functions written for the purpose of the assignment. I have used Django to serve up the code with a visual GUI. The template for the Django server has been adapted based on  the tutorial [How to Render Django Form Manually](https://simpleisbetterthancomplex.com/article/2017/08/19/how-to-render-django-form-manually.html)

Requirements: SRA Toolkit, BWA, Samtools, the whole sacCer3 genome in FASTA format, Django

Key files:
 - [script_cbb520.py](https://github.com/DevangThakkar/Calling_SNPs_Indels/blob/master/script_cbb50.py): Contains all the actual code including shell commands
 - [wg.fa](https://github.com/DevangThakkar/Calling_SNPs_Indels/blob/master/wg.fa): Reference genome for sacCer3
 - [site/manage.py](https://github.com/DevangThakkar/Calling_SNPs_Indels/blob/master/site/manage.py): Settings for Django
 - [site/mysite/forms.py](https://github.com/DevangThakkar/Calling_SNPs_Indels/blob/master/site/mysite/forms.py): Logic behind the visual interface, calls functions from script_cbb520

Questions may be directed to firstname dot lastname at duke dot edu
