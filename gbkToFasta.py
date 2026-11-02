from Bio import SeqIO
gbk_filename = "data/inputs/NM_000513.gbk"
faa_filename = "data/outputs/NM_000513.fas"
input_handle  = open(gbk_filename, "r")
output_handle = open(faa_filename, "w")

i =0;
for seq_record in SeqIO.parse(input_handle, "genbank") :
    print ("Dealing with GenBank record %s" % seq_record.id)
    for seq_feature in seq_record.features :
    # print ("feat %d",i)
        if seq_feature.type=="CDS" :
            assert len(seq_feature.qualifiers['translation'])==1
            output_handle.write(">%s: %s\n%s\n" % (
                   seq_record.name,
                   seq_record.description,
                   seq_feature.qualifiers['translation'][0]))
        # i+=1;
output_handle.close()
input_handle.close()