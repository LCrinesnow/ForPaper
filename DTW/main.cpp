//
//  main.cpp
//  DTW
//
//  Created by Rinesnow on 15/5/6.
//  Copyright (c) 2015年 Rinesnow. All rights reserved.
//

/***********************************************************************/
/************************* DISCLAIMER **********************************/
/***********************************************************************/
/** This UCR Suite software is copyright protected © 2012 by          **/
/** Thanawin Rakthanmanon, Bilson Campana, Abdullah Mueen,            **/
/** Gustavo Batista and Eamonn Keogh.                                 **/
/**                                                                   **/
/** Unless stated otherwise, all software is provided free of charge. **/
/** As well, all software is provided on an "as is" basis without     **/
/** warranty of any kind, express or implied. Under no circumstances  **/
/** and under no legal theory, whether in tort, contract,or otherwise,**/
/** shall Thanawin Rakthanmanon, Bilson Campana, Abdullah Mueen,      **/
/** Gustavo Batista, or Eamonn Keogh be liable to you or to any other **/
/** person for any indirect, special, incidental, or consequential    **/
/** damages of any character including, without limitation, damages   **/
/** for loss of goodwill, work stoppage, computer failure or          **/
/** malfunction, or for any and all other damages or losses.          **/
/**                                                                   **/
/** If you do not agree with these terms, then you you are advised to **/
/** not use this software.                                            **/
/***********************************************************************/
/***********************************************************************/


#include <stdio.h>
#include <stdlib.h>
//#include <math.h>
#include <cmath>
#include <time.h>
#include <iostream>
#include <string>
#include <fstream>
#include<algorithm>
#include <vector>

#define min(x,y) ((x)<(y)?(x):(y))
#define max(x,y) ((x)>(y)?(x):(y))
#define dist(x,y) ((x-y)*(x-y))

#define INF 1e20       //Pseudo Infitinte number for this code

using namespace std;
/*------------------------------
 LC SAX Discord Structs&Methods
 -------------------------------*/
//LC
typedef struct SAX
{
    char word[3];
    int index;
    int count;
} SAX;

const int sonnum=3,base='1';
typedef struct Trie
{
    int num;//记录多少单词途径该节点，即多少单词拥有以该节点为末尾的前缀
    bool terminal;//若terminal==true，该节点没有后续节点
    int count;//记录单词的出现次数，此节点即一个完整单词的末尾字母
    struct Trie *son[sonnum];//后续节点
    vector<int> Index;//索引数组
    
} Trie;

typedef struct DistLoc
{
    double distance;
    int loc;
    
} DL;

/*********************************
 创建一个新节点
 *********************************/
Trie *NewTrie()
{
    Trie *temp=new Trie;
    temp->num=1;
    temp->terminal=false;
    temp->count=0;
    for(int i=0;i<sonnum;++i)temp->son[i]=NULL;
    return temp;
}
/*********************************
 插入一个新词到字典树
 pnt:树根
 s  :新词
 len:新词长度
 *********************************/
void Insert(Trie *pnt,char *s,int len,int position, SAX sax[])
{
    Trie *temp=pnt;
    for(int i=0;i<len;++i)
    {
        if(temp->son[s[i]-base]==NULL){
            temp->son[s[i]-base]=NewTrie();
        }
        else {temp->son[s[i]-base]->num++;temp->terminal=false;}
        temp=temp->son[s[i]-base];
    }
    temp->terminal=true;
    temp->count++;
    temp->Index.push_back(position);
    
    for (int j=0;j<=temp->count-1;j++) {//更新countarray的值，让countarray 做到第0个112和第7个112都是4
        sax[temp->Index.at(j)].count=temp->count;
    }
    
    
}
/*********************************
 删除整棵树
 pnt:树根
 *********************************/
void Delete(Trie *pnt)
{
    if(pnt!=NULL)
    {
        for(int i=0;i<sonnum;++i)if(pnt->son[i]!=NULL)Delete(pnt->son[i]);
        delete pnt;
        pnt=NULL;
    }
}
/*********************************
 查找单词在字典树中的末尾节点
 pnt:树根
 s  :单词
 len:单词长度
 *********************************/
Trie *Find(Trie *pnt,char *s,int len)
{
    Trie *temp=pnt;
    for(int i=0;i<len;++i)
        if(temp->son[s[i]-base]!=NULL){
            temp=temp->son[s[i]-base];
            cout<<"count值"<<temp->count<<endl;
            cout<<"num值"<<temp->num<<endl;
            
        }
        else return NULL;
    //每个被查找的字符串对应的Index们
    for (int j=0;j<=temp->count-1;j++) {
        cout<<"Index:"<<temp->Index.at(j)<<endl;
    }
    return temp;
}

/********************************
 寻找count最小的p集合（SAXword）
 并以此寻找q集合（Index）
 *******************************/
int Outter(int countarray[],int actualength)
{
    
    return 0;
}
/********************************
 qsort的compare函数
 *******************************/

int compareSAX(const void *a, const void* b)
{   SAX* x = (SAX*)a;
    SAX* y = (SAX*)b;
    return x->count - y->count;   // high to low
}
/*------------------------------
 LC SAX Discord Structs&Methods
 -------------------------------*/


/// Data structure for sorting the query
typedef struct Index
{   double value;
    int    index;
} Index;

/// Data structure (circular array) for finding minimum and maximum for LB_Keogh envolop
struct deque
{   int *dq;
    int size,capacity;
    int f,r;
};


/// Sorting function for the query, sort by abs(z_norm(q[i])) from high to low
int comp(const void *a, const void* b)
{   Index* x = (Index*)a;
    Index* y = (Index*)b;
    return abs(y->value) - abs(x->value);   // high to low
}
int Qscomp(const void *a, const void* b)
{
    return ( *(int*)a - *(int*)b );
}


/// Initial the queue at the begining step of envelop calculation
void init(deque *d, int capacity)
{
    d->capacity = capacity;
    d->size = 0;
    d->dq = (int *) malloc(sizeof(int)*d->capacity);
    d->f = 0;
    d->r = d->capacity-1;
}

/// Destroy the queue
void destroy(deque *d)
{
    free(d->dq);
}

/// Insert to the queue at the back
void push_back(struct deque *d, int v)
{
    d->dq[d->r] = v;
    d->r--;
    if (d->r < 0)
        d->r = d->capacity-1;
    d->size++;
}

/// Delete the current (front) element from queue
void pop_front(struct deque *d)
{
    d->f--;
    if (d->f < 0)
        d->f = d->capacity-1;
    d->size--;
}

/// Delete the last element from queue
void pop_back(struct deque *d)
{
    d->r = (d->r+1)%d->capacity;
    d->size--;
}

/// Get the value at the current position of the circular queue
int front(struct deque *d)
{
    int aux = d->f - 1;
    
    if (aux < 0)
        aux = d->capacity-1;
    return d->dq[aux];
}

/// Get the value at the last position of the circular queueint back(struct deque *d)
int back(struct deque *d)
{
    int aux = (d->r+1)%d->capacity;
    return d->dq[aux];
}

/// Check whether or not the queue is empty
int empty(struct deque *d)
{
    return d->size == 0;
}

/// Finding the envelop of min and max value for LB_Keogh
/// Implementation idea is intoruduced by Danial Lemire in his paper
/// "Faster Retrieval with a Two-Pass Dynamic-Time-Warping Lower Bound", Pattern Recognition 42(9), 2009.
void lower_upper_lemire(double *t, int len, int r, double *l, double *u)
{
    struct deque du, dl;
    
    init(&du, 2*r+2);
    init(&dl, 2*r+2);
    
    push_back(&du, 0);
    push_back(&dl, 0);
    
    for (int i = 1; i < len; i++)
    {
        if (i > r)
        {
            u[i-r-1] = t[front(&du)];
            l[i-r-1] = t[front(&dl)];
        }
        if (t[i] > t[i-1])
        {
            pop_back(&du);
            while (!empty(&du) && t[i] > t[back(&du)])
                pop_back(&du);
        }
        else
        {
            pop_back(&dl);
            while (!empty(&dl) && t[i] < t[back(&dl)])
                pop_back(&dl);
        }
        push_back(&du, i);
        push_back(&dl, i);
        if (i == 2 * r + 1 + front(&du))
            pop_front(&du);
        else if (i == 2 * r + 1 + front(&dl))
            pop_front(&dl);
    }
    for (int i = len; i < len+r+1; i++)
    {
        u[i-r-1] = t[front(&du)];
        l[i-r-1] = t[front(&dl)];
        if (i-front(&du) >= 2 * r + 1)
            pop_front(&du);
        if (i-front(&dl) >= 2 * r + 1)
            pop_front(&dl);
    }
    destroy(&du);
    destroy(&dl);
}

/// Calculate quick lower bound
/// Usually, LB_Kim take time O(m) for finding top,bottom,fist and last.
/// However, because of z-normalization the top and bottom cannot give siginifant benefits.
/// And using the first and last points can be computed in constant time.
/// The prunning power of LB_Kim is non-trivial, especially when the query is not long, say in length 128.
double lb_kim_hierarchy(double *t, double *q, int j, int len, double mean, double std, double bsf = INF)
{
    /// 1 point at front and back
    double d, lb;
    double x0 = (t[j] - mean) / std;
    double y0 = (t[(len-1+j)] - mean) / std;
    lb = dist(x0,q[0]) + dist(y0,q[len-1]);
    if (lb >= bsf)   return lb;
    
    /// 2 points at front
    double x1 = (t[(j+1)] - mean) / std;
    d = min(dist(x1,q[0]), dist(x0,q[1]));
    d = min(d, dist(x1,q[1]));
    lb += d;
    if (lb >= bsf)   return lb;
    
    /// 2 points at back
    double y1 = (t[(len-2+j)] - mean) / std;
    d = min(dist(y1,q[len-1]), dist(y0, q[len-2]) );
    d = min(d, dist(y1,q[len-2]));
    lb += d;
    if (lb >= bsf)   return lb;
    
    /// 3 points at front
    double x2 = (t[(j+2)] - mean) / std;
    d = min(dist(x0,q[2]), dist(x1, q[2]));
    d = min(d, dist(x2,q[2]));
    d = min(d, dist(x2,q[1]));
    d = min(d, dist(x2,q[0]));
    lb += d;
    if (lb >= bsf)   return lb;
    
    /// 3 points at back
    double y2 = (t[(len-3+j)] - mean) / std;
    d = min(dist(y0,q[len-3]), dist(y1, q[len-3]));
    d = min(d, dist(y2,q[len-3]));
    d = min(d, dist(y2,q[len-2]));
    d = min(d, dist(y2,q[len-1]));
    lb += d;
    
    return lb;
}

/// LB_Keogh 1: Create Envelop for the query
/// Note that because the query is known, envelop can be created once at the begenining.
///
/// Variable Explanation,
/// order : sorted indices for the query.
/// uo, lo: upper and lower envelops for the query, which already sorted.
/// t     : a circular array keeping the current data.
/// j     : index of the starting location in t
/// cb    : (output) current bound at each position. It will be used later for early abandoning in DTW.
double lb_keogh_cumulative(int* order, double *t, double *uo, double *lo, double *cb, int j, int len, double mean, double std, double best_so_far = INF)
{
    double lb = 0;
    double x, d;
    
    for (int i = 0; i < len && lb < best_so_far; i++)
    {
        x = (t[(order[i]+j)] - mean) / std;
        d = 0;
        if (x > uo[i])
            d = dist(x,uo[i]);
        else if(x < lo[i])
            d = dist(x,lo[i]);
        lb += d;
        cb[order[i]] = d;
    }
    return lb;
}

/// LB_Keogh 2: Create Envelop for the data
/// Note that the envelops have been created (in main function) when each data point has been read.
///
/// Variable Explanation,
/// tz: Z-normalized data
/// qo: sorted query
/// cb: (output) current bound at each position. Used later for early abandoning in DTW.
/// l,u: lower and upper envelop of the current data
double lb_keogh_data_cumulative(int* order, double *tz, double *qo, double *cb, double *l, double *u, int len, double mean, double std, double best_so_far = INF)
{
    double lb = 0;
    double uu,ll,d;
    
    for (int i = 0; i < len && lb < best_so_far; i++)
    {
        uu = (u[order[i]]-mean)/std;
        ll = (l[order[i]]-mean)/std;
        d = 0;
        if (qo[i] > uu)
            d = dist(qo[i], uu);
        else
        {   if(qo[i] < ll)
            d = dist(qo[i], ll);
        }
        lb += d;
        cb[order[i]] = d;
    }
    //    cout<<"检测lb_k2:"<<lb<<endl;
    
    return lb;
}

/// Calculate Dynamic Time Wrapping distance
/// A,B: data and query, respectively
/// cb : cummulative bound used for early abandoning
/// r  : size of Sakoe-Chiba warpping band
double dtw(double* A, double* B, double *cb, int m, int r, double bsf = INF)
{
    
    double *cost;
    double *cost_prev;
    double *cost_tmp;
    int i,j,k;
    double x,y,z,min_cost;
    
    /// Instead of using matrix of size O(m^2) or O(mr), we will reuse two array of size O(r).
    cost = (double*)malloc(sizeof(double)*(2*r+1));
    for(k=0; k<2*r+1; k++)    cost[k]=INF;
    
    cost_prev = (double*)malloc(sizeof(double)*(2*r+1));
    for(k=0; k<2*r+1; k++)    cost_prev[k]=INF;
    
    for (i=0; i<m; i++)
    {
        k = max(0,r-i);
        min_cost = INF;
        
        for(j=max(0,i-r); j<=min(m-1,i+r); j++, k++)
        {
            /// Initialize all row and column
            if ((i==0)&&(j==0))
            {
                cost[k]=dist(A[0],B[0]);
                min_cost = cost[k];
                continue;
            }
            
            if ((j-1<0)||(k-1<0))     y = INF;
            else                      y = cost[k-1];
            if ((i-1<0)||(k+1>2*r))   x = INF;
            else                      x = cost_prev[k+1];
            if ((i-1<0)||(j-1<0))     z = INF;
            else                      z = cost_prev[k];
            
            /// Classic DTW calculation
            cost[k] = min( min( x, y) , z) + dist(A[i],B[j]);
            
            /// Find minimum cost in row for early abandoning (possibly to use column instead of row).
            if (cost[k] < min_cost)
            {   min_cost = cost[k];
            }
        }
        
        /// We can abandon early if the current cummulative distace with lower bound together are larger than bsf
        if (i+r < m-1 && min_cost + cb[i+r+1] >= bsf)
        {   free(cost);
            free(cost_prev);
            return min_cost + cb[i+r+1];
        }
        
        /// Move current array to previous array.
        cost_tmp = cost;
        cost = cost_prev;
        cost_prev = cost_tmp;
    }
    k--;
    
    /// the DTW distance is in the last cell in the matrix of size O(m^2) or at the middle of our array.
    double final_dtw = cost_prev[k];
    free(cost);
    free(cost_prev);
    return final_dtw;
}

/// Print function for debugging
void printArray(double *x, int len)
{   for(int i=0; i<len; i++)
    printf(" %6.2lf",x[i]);
    printf("\n");
}

/// If expected error happens, teminated the program.
void error(int id)
{
    if(id==1)
        printf("ERROR : Memory can't be allocated!!!\n\n");
    else if ( id == 2 )
        printf("ERROR : File not Found!!!\n\n");
    else if ( id == 3 )
        printf("ERROR : Can't create Output File!!!\n\n");
    else if ( id == 4 )
    {
        printf("ERROR : Invalid Number of Arguments!!!\n");
        printf("Command Usage:  UCR_DTW.exe  data-file  query-file   m   R\n\n");
        printf("For example  :  UCR_DTW.exe  data.txt   query.txt   128  0.05\n");
    }
    exit(1);
}


int loc2it(int orign,int m,int EPOCH)
{
    //    return (orign-i+m-1)/EPOCH+m-1;
    return (orign-m+1)/(EPOCH-m+1);
}
int loc2i(int orign,int it,int m,int EPOCH)

{
    return orign-1-(it*(EPOCH-m+1))+m;//减了个m
}






















/// Main Function
int main(  int argc , char *argv[] )
{
    
    int array_size = 1024; // 需要随着不同的datapoint长度随时调整，比如ECG数据Datapoint时3571被11整除以后《1024，可以用define the size of character array
    char * array = new char[array_size]; // allocating an array of 1kb
    int position = 0; //this will be used incremently to fill characters in the array
    
    ifstream fin("/Users/rinesnow/Github/DTW/DTW/FirstSymbol.txt");
    //opening an input stream for file test.txt
    /*checking whether file could be opened or not. If file does not exist or don't have read permissions, file
     stream could not be opened.*/
    if(fin.is_open())
    {
        //file opened successfully so we are here
        cout << "File Opened successfully!!!. Reading data from file into array" << endl;
        //this loop run until end of file (eof) does not occur
        while(!fin.eof() && position < array_size)
        {
            
            fin.get(array[position]); //reading one character from file to array
            position++;
            
        }
        
        array[position-1] = '\0'; //placing character array terminating character
        
        cout << "Displaying Array..." << endl << endl;
        //this loop display all the charaters in array till \0
        for(int ii = 0; array[ii] != '\0'; ii++)
        {
            cout << array[ii];
        }
    }
    else //file could not be opened
    {
        cout << "File could not be opened." << endl;
    }
    cout<<position<<endl;
    //以每三个char为一行
    //   string *wordarray =new string[array_size/3];
    int jj=0;
    int iw=0;//write i
    SAX *sax= new SAX[array_size/3];
    
    while (array[iw]!='\0') {
        sax[jj].word[0]=array[iw];
        sax[jj].word[1]=array[iw+1];
        sax[jj].word[2]=array[iw+2];
        sax[jj].index=jj+1;
        cout<<"排序qian"<<sax[jj].index<<" "<<sax[jj].word[0]<<sax[jj].word[1]<<sax[jj].word[2]<<endl;
        jj++;
        iw+=3;
    }
    cout<<jj<<endl;
    int actualen=jj-1;//wordarray中实际数组的大小
    cout<<actualen<<endl;
    sax[jj].word[0]='E';
    
    /*********************************
     二维数组编写结束
     下面将二维数组录入Trie中
     *********************************/
    Trie *temp;
    Trie *saxtemp;
    temp=NewTrie();
    saxtemp=NewTrie();
    int kk=0;
    while (sax[kk].word[0]!='E') {
        Insert(saxtemp, sax[kk].word, 3, kk,sax);
        kk++;
    }
    /********************************
     sax根据count 由小到大排序
     确定了outter的顺序count的p集合（sax）
     *******************************/
    qsort(sax,actualen,sizeof(SAX),compareSAX);
    int n;
    for (n=0;n<actualen;n++)
        cout<<n<<"排序后"<<sax[n].index<<" "<<sax[n].word[0]<<sax[n].word[1]<<sax[n].word[2]<<" "<<sax[n].count<<endl;
    /********************************
     根据p找q
     开始Discord Detection的代码
     *******************************/
    //  Trie *t;
    //  t=Find(temp,word[0], 3);
    float BSFD=0;
    int compressionrate=16;
    float NND=99999;
    int BSFL;
    for (int P=0; P<actualen;P++) {
        NND=99999;
    }
    /**********
     test
     *********/
    
    
    
    
    cout<<"hello"<<endl;
    //    return 0;
    
    
    FILE *fp,*fp1;            /// data file pointer
    //    FILE *qp;            /// query file pointer
    double bsf;          /// best-so-far
    double *t, *q;       /// data array and query array
    int *order;          ///new order of the query
    double *u, *l, *qo, *uo, *lo,*tz,*cb, *cb1, *cb2,*u_d, *l_d;
    
    
    
    double d;
    long long i,j;
    double ex,ex2,mean,std;
    int m=-1, r=-1;
    long long loc = 0;
    double t1,t2;
    int kim = 0,keogh = 0, keogh2 = 0;
    double dist=0, lb_kim=0, lb_k=0, lb_k2=0;
    double *buffer2, *u_buff, *l_buff,*buffer1;
    Index *Q_tmp;
    
    string dataPath, queryPath;
    
    /// For every EPOCH points, all cummulative values, such as ex (sum), ex2 (sum square), will be restarted for reducing the floating point error.
    int EPOCH = 1000;
    
    
    
    //    cout << "请输入query(m)的大小" <<endl;
    //    cin >> m;
    float R;
    //    cout << "请输入R值"<<endl;
    //    cin >> R;
    R=0.05;
    m=240;// here's m
    if (R<=1)
        r = floor(R*m);
    else
        r = floor(R);
    
    
    
    /// start the clock
    t1 = clock();
    
    
    /// malloc everything here
    //    q = (double *)malloc(sizeof(double)*m);
    //    if( q == NULL )
    //        error(1);
    qo = (double *)malloc(sizeof(double)*m);
    if( qo == NULL )
        error(1);
    uo = (double *)malloc(sizeof(double)*m);
    if( uo == NULL )
        error(1);
    lo = (double *)malloc(sizeof(double)*m);
    if( lo == NULL )
        error(1);
    
    order = (int *)malloc(sizeof(int)*m);
    if( order == NULL )
        error(1);
    //
    //    Q_tmp = (Index *)malloc(sizeof(Index)*m);
    //    if( Q_tmp == NULL )
    //        error(1);
    
    u = (double *)malloc(sizeof(double)*m);
    if( u == NULL )
        error(1);
    
    l = (double *)malloc(sizeof(double)*m);
    if( l == NULL )
        error(1);
    
    cb = (double *)malloc(sizeof(double)*m);
    if( cb == NULL )
        error(1);
    
    cb1 = (double *)malloc(sizeof(double)*m);
    if( cb1 == NULL )
        error(1);
    
    cb2 = (double *)malloc(sizeof(double)*m);
    if( cb2 == NULL )
        error(1);
    
    u_d = (double *)malloc(sizeof(double)*m);
    if( u == NULL )
        error(1);
    
    l_d = (double *)malloc(sizeof(double)*m);
    if( l == NULL )
        error(1);
    
    t = (double *)malloc(sizeof(double)*m*2);
    if( t == NULL )
        error(1);
    
    tz = (double *)malloc(sizeof(double)*m);
    if( tz == NULL )
        error(1);
    
    buffer2 = (double *)malloc(sizeof(double)*EPOCH);
    if( buffer2 == NULL )
        error(1);
    
    buffer1 = (double *)malloc(sizeof(double)*m);
    if( buffer1 == NULL )
        error(1);
    
    u_buff = (double *)malloc(sizeof(double)*EPOCH);
    if( u_buff == NULL )
        error(1);
    
    l_buff = (double *)malloc(sizeof(double)*EPOCH);
    if( l_buff == NULL )
        error(1);
    q = (double *)malloc(sizeof(double)*m);
    if( q == NULL )
        error(1);
    int frontbound;//各个it的前界
    int backbound;//各个it的后界
    bool ifoverwhelm=false;
    int exceptionbackbound;//暂存的特殊it的后界
    
    long int orign;  //数据要时长的话需要扩大成longlongint
    
    DL disloc[actualen];
    
    cout<<"q+orign-1:"<<orign-1<<endl;
    
    for (int outter=0; outter<actualen; outter++) {
        
        fp = fopen("/Users/rinesnow/Github/DTW/DTW/First.txt","r");//若想变成输入式参考上边代码
        if( fp == NULL )
            error(2);

        int countifall=0;
        /// Read query file
        bsf = INF;
        i = 0;
        j = 0;
        ex = ex2 = 0;
        
        
        orign=sax[outter].index*compressionrate*3;
        i=0;
        int w=0;
        //如果文件太大就无法通过将所有数据都放进内存而必须跟buffer2一样
        while(fscanf(fp,"%lf",&d) != EOF)
        {
//            q[i]=d;
//            if(orign+m>)
            if ((i>orign-2)&&(i<orign+m)){
                buffer1[w]= d;
//              cout<<d<<endl;
                w++;
            }
            i++;
        }
    
 
//        memmove(buffer1,q+orign-1,m*sizeof(double));
                cout<<"q+orign-1:"<<orign-1<<endl;
    cout<<"原始序列"<<orign<<" "<<buffer1[0]<<buffer1[1]<<buffer1[2]<<buffer1[77]<<buffer1[78]<<buffer1[79]<<endl;
    
        //        free(q);
        fclose(fp);
        //如果文件太大就无法通过将所有数据都放进内存而必须跟buffer2一样
        i=0;
        while(i<m)
        {
            d = buffer1[i];
                        cout<<i<<" ";
                        cout<<buffer1[i]<<endl;
            
            ex += d;
            ex2 += d*d;
            i++;
        }
        //        cout<<"ex:"<<ex<<"ex2"<<ex2<<endl;
        /// Do z-normalize the query, keep in same array, q
        mean = ex/m;
        std = ex2/m;
        std = sqrt(std-mean*mean);
        for( i = 0 ; i < m ; i++ )
            buffer1[i] = (buffer1[i] - mean)/std;
        
        
        /// Create envelop of the query: lower envelop, l, and upper envelop, u
        lower_upper_lemire(buffer1, m, r, l, u);
        
        Q_tmp = (Index *)malloc(sizeof(Index)*m);
        if( Q_tmp == NULL )
            error(1);
        
        /// Sort the query one time by abs(z-norm(q[i]))
        for( i = 0; i<m; i++)
        {
            Q_tmp[i].value = buffer1[i];
            Q_tmp[i].index = i;
        }
        qsort(Q_tmp, m, sizeof(Index),comp);
        
        /// also create another arrays for keeping sorted envelop
        for( i=0; i<m; i++)
        {   int o = Q_tmp[i].index;
            order[i] = o;
            qo[i] = buffer1[o];
            uo[i] = u[o];
            lo[i] = l[o];
        }
        free(Q_tmp);
        
        /// Initial the cummulative lower bound
        for( i=0; i<m; i++)
        {   cb[i]=0;
            cb1[i]=0;
            cb2[i]=0;
        }
        
        i = 0;          /// current index of the data in current chunk of size EPOCH
        j = 0;          /// the starting index of the data in the circular array, t
        ex = ex2 = 0;
        bool done = false;
        int it=0, ep=0, k=0;
        long long I;    /// the starting index of the data in current chunk of size EPOCH
        
        fp1 = fopen("/Users/rinesnow/Github/DTW/DTW/First.txt","r");//若想变成输入式参考上边代码
        if( fp1 == NULL )
            error(2);
        
        while(!done)
        {
            
            /// Read first m-1 points
            ep=0;
            if (it==0)
            {
                
                for(k=0; k<m-1; k++){
                    
                    if (fscanf(fp1,"%lf",&d)!= EOF){
                        buffer2[k] = d;

                    }
                }
                //                cout<<"序列完毕"<<buffer2[k-1]<<buffer2[k]<<endl;
            }
            else
            {
                for(k=0; k<m-1; k++)
                    buffer2[k] = buffer2[EPOCH-m+1+k];
            }
            
            // Read buffer of size EPOCH or when all data has been read.
            ep=m-1;//1000个1000个读
            while(ep<EPOCH)
            {
                if (fscanf(fp1,"%lf",&d) == EOF)
                    break;
                buffer2[ep] = d;
                ep++;
            }
            //            cout<<"EPOCH完毕"<<buffer2[ep-1]<<buffer2[m-1]<<endl;
            
            // Data are read in chunk of size EPOCH.
            // When there is nothing to read, the loop is end.
            if (ep<=m-1)
            {
                done = true;
            } else {
                
                lower_upper_lemire(buffer2, ep, r, l_buff, u_buff);
                
                // Just for printing a dot for approximate a million point. Not much accurate.
                if (it%(1000000/(EPOCH-m+1))==0)
                    fprintf(stderr,".");
                
                /// Do main task here..
                //                cout<<"轮数："<<it<<endl;
                

                ex=0; ex2=0;
                for(i=0;i<ep;i++)//ep= 整个datafile的长度<EPOCH
                {
                    frontbound=loc2i(orign,it,m,EPOCH)-m-m;
                    backbound=loc2i(orign,it,m,EPOCH)+m-m;
                    
//                    cout<<"backbound"<<backbound<<endl;

                    /// A bunch of data has been read and pick one of them at a time to use
                    d = buffer2[i];
//                    cout<<"it:"<<it;
//                    cout<<" i1:"<<loc2i(orign,it,m,EPOCH)-m-m<<" i2:"<<loc2i(orign,it,m,EPOCH)+m-m;
//                    cout<<" loc:"<<(it)*(EPOCH-m+1)+i-m+1;
//                    cout<<" d: "<<d<<endl;
                   
                        if((it==loc2it(orign, m, EPOCH))&&(i>=frontbound)&&(i<=backbound))//LC实现｜p－q｜>=n
                        {
                            if(backbound-ep>0)
                                ifoverwhelm=true;
                            exceptionbackbound=backbound;
//                            if(countifall<=2*m){
                                countifall++;
                                // 用1021-999+1=21   160-139=21
//                                cout<<"跳跳it:"<<it;
//                                cout<<" i1:"<<frontbound<<" i2:"<<backbound;
//                                cout<<" loc:"<<(it)*(EPOCH-m+1)+i-m+1;
//                                cout<<" d: "<<d;
//                                cout<<" count:"<<countifall<<endl;
                                continue;
//                            }
                        }
                    
                        if (it==(loc2it(orign,m,EPOCH)+1)&&(ifoverwhelm==true)&&(countifall<2*m))
                        {
                            if ((i>=m-1)&&(i<=exceptionbackbound-EPOCH+m-1)){
                                countifall++;
                                // 用1021-999+1=21   160-139=21
//                                cout<<"1000后it:"<<it;
//                                cout<<" i1:"<<frontbound<<" i2:"<<backbound;
//                                cout<<" loc:"<<(it)*(EPOCH-m+1)+i-m+1;
//                                cout<<" d: "<<d;
//                                cout<<" count:"<<countifall<<endl;
                                continue;
                            }
                        }

                    
                    
                    
                    
                    
                    /// Calcualte sum and sum square
                    ex += d;
                    ex2 += d*d;
                    
                    /// t is a circular array for keeping current data
                    t[i%m] = d;
                    
                    /// Double the size for avoiding using modulo "%" operator
                    t[(i%m)+m] = d;
                    
                    /// Start the task when there are more than m-1 points in the current chunk
                    

                    if( i >= m-1 )//++
                    {
//                        cout<<i<<endl;

                        mean = ex/m;
                        std = ex2/m;
                        std = sqrt(std-mean*mean);//LC改，因为数据有负的所以不改成绝对值，不能进行下去了
                        /// compute the start location of the data in the current circular array, t
                        j = (i+1)%m;
                        /// the start location of the data in the current chunk
                        I = i-(m-1);
                        /// Use a constant lower bound to prune the obvious subsequence
                        lb_kim = lb_kim_hierarchy(t, buffer1, j, m, mean, std, bsf);
//                        cout<<"lbkim:"<<lb_kim<<endl;
                        if (lb_kim < bsf)
                        {
                            /// Use a linear time lower bound to prune; z_normalization of t will be computed on the fly.
                            /// uo, lo are envelop of the query.
                            lb_k = lb_keogh_cumulative(order, t, uo, lo, cb1, j, m, mean, std, bsf);
//                            cout<<"lbkim:"<<lb_k<<endl;
                            if (lb_k < bsf)
                            {
                                /// Take another linear time to compute z_normalization of t.
                                /// Note that for better optimization, this can merge to the previous function.
                                for(k=0;k<m;k++)
                                {
                                    tz[k] = (t[(k+j)] - mean)/std;
//                                    cout<<tz[k];
                                }
                                
                                /// Use another lb_keogh to prune
                                /// qo is the sorted query. tz is unsorted z_normalized data.
                                /// l_buff, u_buff are big envelop for all data in this chunk
                                lb_k2 = lb_keogh_data_cumulative(order, tz, qo, cb2, l_buff+I, u_buff+I, m, mean, std, bsf);
                                if (lb_k2 < bsf)
                                {
                                    /// Choose better lower bound between lb_keogh and lb_keogh2 to be used in early abandoning DTW
                                    /// Note that cb and cb2 will be cumulative summed here.
                                    if (lb_k > lb_k2)
                                    {
                                        cb[m-1]=cb1[m-1];
                                        for(k=m-2; k>=0; k--)
                                            cb[k] = cb[k+1]+cb1[k];
                                    }
                                    else
                                    {
                                        cb[m-1]=cb2[m-1];
                                        for(k=m-2; k>=0; k--)
                                            cb[k] = cb[k+1]+cb2[k];
                                    }
                                    
                                    /// Compute DTW and early abandoning if possible
                                    dist = dtw(tz, buffer1, cb, m, r, bsf);
//                                    cout<<"dist:"<<dist<<endl;

                                    if( dist < bsf )
                                    {   /// Update bsf
                                        /// loc is the real starting location of the nearest neighbor in the file
                                        bsf = dist;
                                        loc = (it)*(EPOCH-m+1) + i-m+1;
//                                                                            cout<<"loc:::"<<loc<<" ";
//                                                                            cout<<dist<<endl;
                                        
                                    }
                                } else
                                    keogh2++;
                            } else
                                keogh++;
                        } else
                            kim++;
                        
                        /// Reduce obsolute points from sum and sum square
                        ex -= t[j];
                        ex2 -= t[j]*t[j];
                    }
                }
                
                /// If the size of last chunk is less then EPOCH, then no more data and terminate.
                if (ep<EPOCH)
                    done=true;
                else
                    it++;
            }
            
            
        }
        i = (it)*(EPOCH-m+1) + ep;
        fclose(fp1);
        disloc[outter].distance=sqrt(bsf);
        disloc[outter].loc=orign;
        
    }

    free(u);
    free(l);
    free(uo);
    free(lo);
    free(qo);
    free(cb);
    free(cb1);
    free(cb2);
    free(tz);
    free(t);
    free(l_d);
    free(u_d);
    //  free(q);
    free(l_buff);
    free(u_buff);
    
    t2 = clock();
    
    for (int test=0; test<actualen; test++) {
                if(disloc[test].distance>1.5)
        cout<<disloc[test].distance<<" "<<disloc[test].loc<<endl;
    }

    
    printf("\n");
    
    
    /// Note that loc and i are long long.
//    cout << "Location : " << loc << endl;
//    cout << "Distance : " << sqrt(bsf) << endl;
//    cout << "Data Scanned : " << i << endl;
//    cout << "Total Execution Time : " << (t2-t1)/CLOCKS_PER_SEC << " sec" << endl;
    
    /// printf is just easier for formating ;)
//    printf("\n");
//    printf("Pruned by LB_Kim    : %6.2f%%\n", ((double) kim / i)*100);
//    printf("Pruned by LB_Keogh  : %6.2f%%\n", ((double) keogh / i)*100);
//    printf("Pruned by LB_Keogh2 : %6.2f%%\n", ((double) keogh2 / i)*100);
//    printf("DTW Calculation     : %6.2f%%\n", 100-(((double)kim+keogh+keogh2)/i*100));
    return 0;
    
}
