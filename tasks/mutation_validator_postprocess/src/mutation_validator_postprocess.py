import sys
import csv

input_maf_fn = sys.argv[1]
output_maf_fn = sys.argv[2]


#assumes no comment rows
ifid = open(input_maf_fn)
dictreader = csv.DictReader(ifid, dialect='excel-tab')

input_fieldnames = dictreader.fieldnames
old_output_fieldnames = []
new_output_fieldnames = []
validation_data_types = []
for fieldname in input_fieldnames:
    
    if fieldname.startswith('discovery_'):
        continue
    old_output_fieldnames.append(fieldname)
    if fieldname.startswith('validation_judgement'):
        fieldname_list = fieldname.split('_')
        data_type = fieldname_list[-1]
        validation_data_types.append(data_type)
        new_fieldname = 'mutval_status_' + data_type
        new_output_fieldnames.append(new_fieldname)

output_fieldnames = old_output_fieldnames + new_output_fieldnames


ofid = open(output_maf_fn,'w')
dictwriter = csv.DictWriter(ofid, fieldnames=output_fieldnames, extrasaction='ignore' ,dialect='excel-tab', lineterminator='\n')
dictwriter.writeheader()

for line in dictreader:
    for data_type in validation_data_types:
        judgement_str = line['validation_judgement_' + data_type]
        try:
            judgement_int = int(judgement_str)
        except:
            judgement_int = -1
        power_str = line['validation_power_' + data_type]
        try:
            power_float = float(power_str)
        except:
            power_float = -1.0
        if judgement_int == -1 or judgement_int == 0:
            judgement_label = 'unvalidated'
        elif judgement_int == 1:
            judgement_label = 'validated'
        elif judgement_int == 2:
            judgement_label = 'germline'
        else:
            raise Exception('unexpected value for judgement: %s %d'%(judgement_str, judgement_int))
        if power_float >= 0.95:
            power_label = 'powered'
        else:
            power_label = 'unpowered'

        label_field = 'mutval_status_' + data_type
        line[label_field] = judgement_label + '_' + power_label
    dictwriter.writerow(line)
