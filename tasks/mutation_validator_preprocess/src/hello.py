import sys

salutation = sys.argv[1]
name_file = sys.argv[2]

infid = open(name_file)
name_contents = infid.read()
infid.close()
name = name_contents.strip()


outfid = open('greeting.txt','w')
outfid.write('%s %s\n'%(salutation,name))
outfid.close()
