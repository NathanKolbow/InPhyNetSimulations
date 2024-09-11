with open("snaq.tab", "r") as f:
    lines = f.readlines()
    n500r1_lines = [line for line in lines if line.split(",")[0] == "500" and (line.split(",")[1] == "1" or line.split(",")[1] == "2")]
    
    with open("snaq_n500r12.tab", "w+") as new_f:
        for line in n500r1_lines:
            new_f.write(line)