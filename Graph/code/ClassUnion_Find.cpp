class union_find
{
public:
	union_find(int _len)
	{
		opn=0;
		len=_len+5;
		data=NULL;
		data=new int[len];
		//if(data==NULL)throw std::runtime_error("ERROR!");
		for(int i=1;i<len;i++)data[i]=i;
	}
	~union_find()
	{
		delete[] data;
	}
	int find(int x)
	{
		//if(x>len)throw std::runtime_error("ERROR!");
		return x==data[x]?x:data[x]=find(data[x]);
	}
	void merge(int x,int y)
	{
		x=find(x);
		y=find(y);
		if(x!=y)
		{
			data[x]=y;
			opn++;
		}
	}
	int opn;  //操作数
private:
	int *data;
	int len;
};

/** 建一个长度为len的并查集： union_find(len);
 *  搜索n在哪个集合内： object.find(n);
 *  合并n和m所在的集合： object.find(n,m);
 *  opn存储了合并了多少个集合（并非运行了多少次merge(int,int)函数）
 * */


