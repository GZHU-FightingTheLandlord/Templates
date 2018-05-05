#define USE_ZKWTREE
#ifdef USE_ZKWTREE
typedef long long int64;
int64 how_to_update(int64 a,int64 b)
{
	return a+b;
}
/** such as max(a,b) min(a,b) ... **/
class zkw_tree
{
public:
	zkw_tree(int _len)
	{
		int t=1;
		while(_len>=t)t*=2;
		len=t*2;
		data=new int64 [len*2+5];
		data[0]=0;
		//if(data==NULL)throw std::runtime_error("Memory Limit Exceeded or some other error.");
	}
	~zkw_tree()
	{
		delete[] data;
	}
	void update(int64 i,int64 x)
	{
		data[i+=len]+=x;
		for(i>>=1;i;i>>=1)data[i]=how_to_update(data[i<<1],data[i<<1|1]);
	}
	int64 query(int s,int t)/// must make sure that s<t
	{
		int64 ans=0; /** there is a point to recompose while you are not update a summation. **/
		for(s+=len-1,t+=len+1;s^t^1;s>>=1,t>>=1)
		{
			if(~s&1)ans=how_to_update(ans,data[s^1]);
			if(t&1)ans=how_to_update(ans,data[t^1]);
		}
		return ans;
	}
private:
	int len;
	int64 *data;
};
#endif // USE_ZKWTREE
/** This code is not to be used yet.
	Welcome to hack it. **/

/** create a zkw_tree: zkw_tree treea(1000);
	下标从1开始。the subscript is from 1.
	更新树的下标为location为value: object.zkw_tree::update(location,value);
	查询区间[from,to]: object.zkw_tree::query(from,to)
**/
