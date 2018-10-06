#include<bits/stdc++.h>
#define loop(i,a,b) for(int i=a;i<b;i++)
#define loopb(i,a,b) for(int i=a;i>=b;i--)
#define loopm(i,a,b,step) for(int i=a;i<b;i+=step)
#define loopbm(i,a,b,step) for(int i=a;i>=b;i-=step)
#define pb(a) push_back(a)
#define mp(a,b) make_pair(a,b)
#define init(arr,val) memset(arr,val,sizeof(arr))
#define INF 1000000007
#define MOD 1000000007
#define BINF 1000000000000000001
#define int long long int
#define double long double

using namespace std;

const int N=2e5+5;
const int LGN=20;
vector<pair<int,int> >g[N];
int dp[LGN][N];
int height[N];
int dist[N];

void dfs(int s,int p=0,int ht=0,int d=0)
{
	dp[0][s]=p;
	loop(i,1,LGN) dp[i][s]=dp[i-1][dp[i-1][s]];
	height[s]=ht;
	dist[s]=d;
	for(auto v:g[s])
	{
		if(v.first!=p)
		{
			dfs(v.first,s,ht+1,d+v.second);
		}
	}
}

//HLD PART

int sub[N];

void dfs1(int s,int p=0)
{
	sub[s]=1;
	for(auto &v:g[s])
	{
		if(v.first!=p)
		{
			dfs1(v.first,s);
			sub[s]+=sub[v.first];
			if(sub[v.first]>sub[g[s][0].first]) swap(v,g[s][0]);
		}
	}
}

int T=0;
int in[N],out[N],rin[N],nxt[N];

void dfs_hld(int s,int p=0)
{
	in[s]=T++;
	rin[in[s]]=s;
	if(s==1) nxt[s]=s;
	for(auto v:g[s])
	{
		if(v.first!=p)
		{
			nxt[v.first]=(v.first==g[s][0].first?nxt[s]:v.first);
			dfs_hld(v.first,s);
		}
	}
}

//LCA

int lca(int u,int v)
{
	if(height[u]<height[v])
	{
		swap(u,v);
	}
	
	int k=height[u]-height[v];
	loop(i,0,LGN) if(k&(1<<i)) u=dp[i][u];
	
	if(u==v) return u;
	
	loopb(i,LGN-1,0)
	{
		if(dp[i][u]!=dp[i][v]) u=dp[i][u],v=dp[i][v];
	}
	return dp[0][u];
}

int length(int u,int v)
{
	//Assuming u is an ancestor of v
	return (dist[v]-dist[u]);
}


/*DYNAMIC CONVEX HULL TRICK*/
/*credits:https://ideone.com/fAjoZu
	   http://codeforces.com/blog/entry/11339*/

/*
 * Dynamic version of data structure
 * to be used in dynamic programming optimisation
 * called "Convex Hull trick"
 * 'Dynamic' means that there is no restriction on adding lines order
 */
 

class ConvexHullDynamic
{
    typedef double coef_t;
    typedef double coord_t;
    typedef double val_t;

    /*
     * Line 'y=a*x+b' represented by 2 coefficients 'a' and 'b'
     * and 'xLeft' which is intersection with previous line in hull(first line has -INF)
     */
private:
    struct Line
    {
        coef_t a, b;
        double xLeft;

        enum Type {line, maxQuery, minQuery} type;
        coord_t val;

        explicit Line(coef_t aa=0, coef_t bb=0) : a(aa), b(bb), xLeft(-INFINITY), type(Type::line), val(0) {}
        val_t valueAt(coord_t x) const { return a*x+b; }
        friend bool areParallel(const Line& l1, const Line& l2) { return l1.a==l2.a; }
        friend double intersectX(const Line& l1, const Line& l2) { return areParallel(l1,l2)?INFINITY:1.0*(l2.b-l1.b)/(l1.a-l2.a); }
        bool operator<(const Line& l2) const
        {
            if (l2.type == line)
                return this->a     < l2.a;
            if (l2.type == maxQuery)
                return this->xLeft < l2.val;
            if (l2.type == minQuery)
                return this->xLeft > l2.val;
        }
    };

private:
    bool            isMax; //whether or not saved envelope is top(search of max value)
    std::set<Line>  hull;  //envelope itself

private:
    /*
     * INFO:        Check position in hull by iterator
     * COMPLEXITY:  O(1)
     */
    bool hasPrev(std::set<Line>::iterator it) { return it!=hull.begin(); }
    bool hasNext(std::set<Line>::iterator it) { return it!=hull.end() && std::next(it)!=hull.end(); }

    /*
     * INFO:        Check whether line l2 is irrelevant
     * NOTE:        Following positioning in hull must be true
     *              l1 is next left to l2
     *              l2 is right between l1 and l3
     *              l3 is next right to l2
     * COMPLEXITY:  O(1)
     */
    bool irrelevant(const Line& l1, const Line& l2, const Line& l3) { return intersectX(l1,l3) <= intersectX(l1,l2); }
    bool irrelevant(std::set<Line>::iterator it)
    {
        return hasPrev(it) && hasNext(it)
               && (    isMax && irrelevant(*std::prev(it), *it, *std::next(it))
                       || !isMax && irrelevant(*std::next(it), *it, *std::prev(it)) );
    }

    /*
     * INFO:        Updates 'xValue' of line pointed by iterator 'it'
     * COMPLEXITY:  O(1)
     */
    std::set<Line>::iterator updateLeftBorder(std::set<Line>::iterator it)
    {
        if (isMax && !hasPrev(it) || !isMax && !hasNext(it))
            return it;

        double val = intersectX(*it, isMax?*std::prev(it):*std::next(it));
        Line buf(*it);
        it = hull.erase(it);
        buf.xLeft = val;
        it = hull.insert(it, buf);
        return it;
    }

public:
    explicit ConvexHullDynamic(bool isMax=false): isMax(isMax) {}

    /*
     * INFO:        Adding line to the envelope
     *              Line is of type 'y=a*x+b' represented by 2 coefficients 'a' and 'b'
     * COMPLEXITY:  Adding N lines(N calls of function) takes O(N*log N) time
     */
    void addLine(coef_t a, coef_t b)
    {
        //find the place where line will be inserted in set
        Line l3 = Line(a, b);
        auto it = hull.lower_bound(l3);

        //if parallel line is already in set, one of them becomes irrelevant
        if (it!=hull.end() && areParallel(*it, l3))
        {
            if (isMax && it->b < b || !isMax && it->b > b)
                it = hull.erase(it);
            else
                return;
        }

        //try to insert
        it = hull.insert(it, l3);
        if (irrelevant(it)) { hull.erase(it); return; }

        //remove lines which became irrelevant after inserting line
        while (hasPrev(it) && irrelevant(std::prev(it))) hull.erase(std::prev(it));
        while (hasNext(it) && irrelevant(std::next(it))) hull.erase(std::next(it));

        //refresh 'xLine'
        it = updateLeftBorder(it);
        if (hasPrev(it))
            updateLeftBorder(std::prev(it));
        if (hasNext(it))
            updateLeftBorder(std::next(it));
    }
    /*
     * INFO:        Query, which returns max/min(depends on hull type - see more info above) value in point with abscissa 'x'
     * COMPLEXITY:  O(log N), N-amount of lines in hull
     */
    val_t getBest(coord_t x) const
    {
        Line q;
        q.val = x;
        q.type = isMax ? Line::Type::maxQuery : Line::Type::minQuery;
        
        //cout<<hull.size()<<endl;
	if(hull.empty()) return BINF;
	
        auto bestLine = hull.lower_bound(q);
        if (isMax) --bestLine;
        return bestLine->valueAt(x);
    }
    
    void clean()
    {
    	hull.clear();
    }
};


//SQUARE ROOT DECOMPOSITION

const int sqn=3;
int belongsto[N];
int block_start[sqn+5];  //[block_start,block_end] in same block
int block_end[sqn+5];
int numblocks=0;
double ans[N];
ConvexHullDynamic cht_ltor[sqn+5];
ConvexHullDynamic cht_rtol[sqn+5];
void decompose()
{
	for(int i=0;i<T;i+=sqn)
	{
		for(int j=i;j<min(T,i+sqn);j++)
		 belongsto[j]=numblocks;
		
		block_start[numblocks]=i;
		block_end[numblocks]=min(i+sqn,T)-1;
		numblocks++;
	}
}

void make_update(int speed,double c,int blkno,int l,int r,int dir)
{
	//Update in block blkno from [l,r]
	//dir=0 ltor means update from u to lca , dir=1 means update from lca to u
	//check if whole block is to be updated
	
	double m=1.0/(double)speed;
	if(l==block_start[blkno] and r==block_end[blkno])
	{
	
		
		if(dir==0)
		{
			cht_ltor[blkno].addLine(m,c);
		}
		else cht_rtol[blkno].addLine(m,c);
	}
	
	else{
		//Partial update. update whole block
		if(dir==0)
		{
			
			//Left to right
			int dd=0;
			while(r>=l)
			{
				//cout<<rin[r]<<" ";
				double newans=c+(double)dd/(double)speed;
				//cout<<fixed<<setprecision(9)<<newans<<endl;
				ans[r]=min(ans[r],newans);
				if(r>l) dd+=length(rin[r-1],rin[r]);
				r--;
			}
		}
		else
		{
			//cout<<l<<" "<<r<<endl;
			int dd=0;
			while(r>=l)
			{
				//cout<<rin[r]<<" ";
				double newans=c-(double)dd/(double)speed;
				//cout<<fixed<<setprecision(9)<<newans<<endl;
				ans[r]=min(ans[r],newans);
				//cout<<fixed<<setprecision(9)<<ans[r]<<endl;
				if(r>l) dd+=length(rin[r-1],rin[r]);
				r--;
			}
		}
		}
	
}


void update_query(int l,int r,int speed,double c,int dir)
{
	int stblk=belongsto[l];
	int enblk=belongsto[r];
	
	
	
        if(dir==0)
        {
        	while(stblk!=enblk)
        	{
        		//cout<<"machao";
        		int beg=block_start[enblk];
        		make_update(speed,c,enblk,beg,r,dir);
        		c+=((double)length(rin[block_end[enblk-1]],rin[r]))/(double)speed;
        		enblk-=1;
        		r=block_end[enblk];
        	}
        	make_update(speed,c,enblk,l,r,dir);
        }
        else
        {
        	//cout<<"YASS\n";
        	while(stblk!=enblk)
        	{
        		//cout<<"YASS\n";
        		int beg=block_start[enblk];
        		make_update(speed,c,enblk,beg,r,dir);
        		c-=((double)length(rin[block_end[enblk-1]],rin[r]))/(double)speed;
        		enblk-=1;
        		r=block_end[enblk];
        	}
        	make_update(speed,c,enblk,l,r,dir);
        }
}


void update(int u,int v,int speed,double c)
{
	int l=lca(u,v);
	int U=u;
	int V=v;
	int L=l;
	//cout<<u<<" "<<v<<" "<<l<<endl;
	//cout<<nxt[8]<<" "<<nxt[l]<<endl;
	double C=c;
	//u to l->dir=0
	while(nxt[u]!=nxt[l])
	{
	//cout<<"PROO HERE\n";
		//cout<<"YASSS  "<<u<<" "<<l<<endl;
		int head=nxt[u];
		update_query(in[head],in[u],speed,C,0);
		int pp=dp[0][head];
		C+=(double)length(pp,u)/(double)speed;
		u=pp;
	}
	//cout<<u<<" "<<"machao  ";
	//cout<<fixed<<setprecision(9)<<C<<endl;
	update_query(in[l],in[u],speed,C,0);
	
	//l to v -> dir=1
	//cout<<length(L,U)<<endl;
	c+=(double)length(L,U)/(double)speed;
	//cout<<fixed<<" "<<setprecision(9)<<c<<endl;
	c+=(double)length(L,V)/(double)speed;
	//cout<<fixed<<" "<<setprecision(9)<<c<<endl;
	while(nxt[v]!=nxt[l])
	{
		//cout<<"YASS\n";
		int head=nxt[v];
		update_query(in[head],in[v],speed,c,1);
		int pp=dp[0][head];
		c-=(double)length(pp,v)/(double)speed;
		v=pp;
	}
	update_query(in[l],in[v],speed,c,1);
}


void fixans(int idx)
{
		
	int blk=belongsto[idx];
	int head=rin[block_start[blk]];
	int tail=rin[block_end[blk]];
	
	int diss=length(rin[idx],tail);
	double ans1=cht_ltor[blk].getBest((double)diss);
	double ans2=cht_rtol[blk].getBest((-1.0*diss));
	//cout<<rin[idx]<<" "<<ans1<<" "<<ans2<<endl;
	ans[idx]=min(ans[idx],min(ans1,ans2));
	
		
}
	

void fuck_everything()
{
	T=0;
	numblocks=0;
	loop(i,0,N)
	{
		g[i].clear();
		height[i]=0;
		dist[i]=0;
		in[i]=out[i]=rin[i]=nxt[i]=0;
		belongsto[i]=0;
		ans[i]=BINF;
		loop(j,0,LGN) dp[j][i]=0;
	}
	loop(i,0,sqn+5)
	{
		block_start[i]=block_end[i]=0;
		cht_ltor[i].clean();
		cht_rtol[i].clean();
	}
}

#undef int
int main()
{
#define int long long int
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    cout.tie(NULL);
    freopen("input.txt","r",stdin);
    int t;
    cin>>t;
    while(t--)
    {
    	fuck_everything();
    	int n;
    	cin>>n;
    	loop(i,0,n-1)
    	{
    		int u,v,w;
    		cin>>u>>v>>w;
    		g[u].push_back({v,w});
    		g[v].push_back({u,w});
    	}
    	
    	dfs(1);
    	dfs1(1);
    	dfs_hld(1);
    	
    	decompose();
    	/*loop(i,1,n+1) cout<<in[i]<<" ";
    	cout<<endl;
    	loop(i,0,T) cout<<rin[i]<<" ";
    	cout<<endl;*/
    	int m;
    	cin>>m;
    	while(m--)
    	{
    		int u,v,con,speed;
    		cin>>u>>v>>con>>speed;
    		
    		double c=1.0*con;
    		double m=1.0/(1.0*speed);
    		//cout<<fixed<<setprecision(9)<<m<<endl;
    		update(u,v,speed,c);
    	}
    	
    	loop(i,0,T) fixans(i);
    	
    	loop(i,1,n+1)
    	{
    		int idx=in[i];
    		double myans=ans[idx];
    		if(myans<BINF) cout<<fixed<<setprecision(9)<<myans<<endl;
    		else cout<<-1<<endl;
    	}
    	
    }
    
    return 0;
}
