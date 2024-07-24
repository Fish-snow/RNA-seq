#open the works
import pandas as pd
up = pd.read_csv("/Users/yuxueqian/Desktop/up.txt",sep='\t')
down=pd.read_csv("/Users/yuxueqian/Desktop/down.txt",sep='\t')
import os
up.to_csv(os.path.join("/Users/yuxueqian/Desktop", "up.csv"))
down.to_csv(os.path.join("/Users/yuxueqian/Desktop", "down.csv"))

#get the target dataframe
up={"Term":up.Term,"Genes": up.Genes}
up=pd.DataFrame(up)
down={"Term":down.Term,"Genes": down.Genes}
down=pd.DataFrame(down)

#split