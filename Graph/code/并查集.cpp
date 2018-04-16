
class union_find
{
public:
	union_find(int _len)
	{
		len=_len;
		data=new int [len];
		init();
	}
	void init()
	{
		for(int i=0;i<len;i++)
		{
			data[i]=i;
		}
	}
	int find(int x)
	{
		return x==data[x]?x:data[x]=find(data[x]);
	}
	void uni(int x,int y)
	{
		x=find(x);
		data[x]=y;
	}
private:
	int len;
	int * data;
};

///Ã÷ÌìÔÙ¸ü
