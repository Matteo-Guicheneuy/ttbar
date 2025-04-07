import sys

def FixHwU(path):
    path2="/home/matteo/mg5amcnlo/"+path+"/FixedOrderAnalysis/"
    file = open(path2+"HwU.f", "r")
    replaced_content = ""
    word=" parameter :: wgts_info_len="
    
    for line in file:
        if(line.find(word)>=0):
            replaced_content+=line.replace("80","50")
        else:
            replaced_content += line 
            
    file.close()
    
    write_file = open(path2+"HwU.f", "w")
    write_file.write(replaced_content)
    write_file.close()


def FixLHEF(path):
    path3="/home/matteo/mg5amcnlo/"+path+"/MCatNLO/include/"
    file = open(path3+"LHEFRead.h", "r")
    replaced_content = ""
    word="const int wgts_info_len_used ="
    
    for line in file:
        if(line.find(word)>=0):
            replaced_content+=line.replace("80","50")
        else:
            replaced_content += line 
            
    file.close()
    
    write_file = open(path3+"LHEFRead.h", "w")
    write_file.write(replaced_content)
    write_file.close()

    

def FixPythia(path):
    path4="/home/matteo/mg5amcnlo/"+path+"/MCatNLO/srcPythia8/"
    file = open(path4+"Pythia83.cc", "r")
    replaced_content = ""
    word="char(*)[50]"
    word2="info[1024][50]"
    
    for line in file:
        if(line.find(word)>=0 or line.find(word2)>=0):
            newline=line.replace("[50]","[80]")
            replaced_content+=newline
        else:
            replaced_content += line 

    file.close()
    
    write_file = open(path4+"Pythia83.cc", "w")
    write_file.write(replaced_content)
    write_file.close()


    
path=sys.argv[1]
#FixHwU(path)
#FixLHEF(path)
FixPythia(path)
