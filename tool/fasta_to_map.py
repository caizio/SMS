# 将fasta格式转成map格式

# 从文件中读取序列id
def readName(path):
    result = []
    with open(path,"r") as file:
        line = file.readline()
        while line:
            words = line.replace("\n","").split(",")
            if len(words) == 2:
                words.append("")
            result.append(words)
            line = file.readline()
    return result

def readFasta(path):
    result = {}
    with open(path,"r") as f:
        name = ""
        seq = ""
        line = f.readline()
        while line:
            if line[0] == ">":
                if seq != "":
                    result[name] = seq
                    seq = ""
                name = line.split("|")[1]
            else:
                seq += line.replace("\n","")
            line = f.readline()
        if seq != "":
            result[name] = seq
    return result

def save_data(path:str,datas:list):
    with open(path,"w") as file:
        for data in datas:
            file.write(data[0] + "," + data[1] + "," + data[2] + "\n")

def fasta_to_map(fasta_path:str,map_path:str,sava_path:str):
    fasta = readFasta(fasta_path)
    map = readName(map_path)
    for line in map:
        if fasta.get(line[1]) != None:
            line[2] = fasta[line[1]]
    save_data(sava_path,map)

if __name__ == "__main__":
    dir = r"D:\code\edit\c++\SMS_new\tool\getSequence\downloads\D. Melanogaster.fasta"
    dir2 = r"D:\code\edit\c++\SMS_new\tool\getSequence\data\D. MelanogasterMap2.txt"
    save = r"D:\code\edit\c++\SMS_new\tool\getSequence\downloads\seq_D. MelanogasterMap2.txt"
    # test = readFasta(dir)
    # print(test)
    # test2 = readName(dir2)
    # print(test2)
    fasta_to_map(dir,dir2,save)