import os
from datetime import datetime
from PyPDF2 import PdfFileWriter, PdfFileReader

# input
input_dir = './writing/'
infile_name = 'tese.pdf'
infile = os.path.join(input_dir, infile_name)

# output
prefix = 'Thesis_Ronaldo_Chaves'
suffix = '.pdf'
output_dir = '.'
outfile_name = prefix + datetime.now().strftime("_%d_%m_%y") + suffix
outfile = os.path.join(output_dir, outfile_name)

# Cut and Trim
input_pdf = PdfFileReader(open(infile, "rb"))
output = PdfFileWriter()

output_file = open(outfile, "wb")
output.addPage(input_pdf.getPage(0))

for i in range(12, input_pdf.getNumPages() - 2):
    output.addPage(input_pdf.getPage(i))

output.write(output_file)
