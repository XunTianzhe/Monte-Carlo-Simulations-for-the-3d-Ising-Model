#include <bits/stdc++.h>
#include <iostream>
#include <random>

using namespace std;

mt19937 gen(std::chrono::system_clock::now().time_since_epoch().count());

int ***canvas;      //整个系统
int ***canvassec;      //测试平衡时间时采用，第二个系统
double Jfactor = 1;
int dx[] = {0, 0, 0, 0, 1, -1};
int dy[] = {0, 0, 1, -1, 0, 0};
int dz[] = {1, -1, 0, 0, 0, 0};
double Tc = 4.511523;

struct p //点结构
{
    int x, y, z;

    p(int x, int y, int z) {
        this->x = x;
        this->y = y;
        this->z = z;
    };
};

//Wolff cluster算法
void clusterflip(int ***&a, double temperature,int SIZE) {

    double pdd; //接受率
    pdd = 1 - exp(-2 * Jfactor / temperature);

    //生成随机数
    uniform_int_distribution<int> u(0, SIZE*SIZE*SIZE-1);
    uniform_real_distribution<double> v(0,1);


    //随机选取一个点作为种子
    int po = u(gen);
    int seedx = po/SIZE/SIZE;
    int seedy = po/SIZE%SIZE;
    int seedz = po%SIZE;
    p seed(seedx, seedy, seedz);
    queue<p> bag; //bag 用来存放种子
    bag.emplace(seed);
    vector<p> cluster; //cluster 用来存放要翻转的点
    cluster.emplace_back(seed);
    bool vis[SIZE][SIZE][SIZE]; //判断一个点是否访问过
    for (int i = 0; i < SIZE; i++)
        for (int j = 0; j < SIZE; j++)
            for (int k = 0; k < SIZE; k++)
                vis[i][j][k] = false;
    vis[seed.x][seed.y][seed.z] = true;

    while (!bag.empty()) {
        p now = bag.front();
        bag.pop();
        for (int i = 0; i < 6; i++) //寻找临域
        {
            int x1 = (now.x + dx[i] + SIZE) % SIZE, y1 = (now.y + dy[i] + SIZE) % SIZE, z1 =
                    (now.z + dz[i] + SIZE) % SIZE;
            double prob = v(gen);
            //如果邻域中的这个点没有被访问过 且 这个点跟种子值同号 同时接受率小于pdd 则将其作为新的种子放入bag中 同时也放入cluster中
            if (!vis[x1][y1][z1] && prob < pdd && a[x1][y1][z1] * a[now.x][now.y][now.z] > 0.1) {
                vis[x1][y1][z1] = true;
                cluster.emplace_back(p(x1, y1, z1));
                bag.emplace(p(x1, y1, z1)); //将这一次的种子丢弃
            }
        }
    }

    for (auto n: cluster)
        a[n.x][n.y][n.z] *= -1;
}

//计算平均磁化强度（把每个点的值相加除以总点数）
double AvrM(int ***&a,int SIZE) {
    int Mag = 0;
    for (int i = 0; i < SIZE; i++)
        for (int j = 0; j < SIZE; j++)
            for (int k = 0; k < SIZE; k++)
                Mag += a[i][j][k];

    return (double) abs(Mag) / (SIZE * SIZE * SIZE);
}

//计算平均能量  即把每个点的能量相加/SIZE^3  每个点的能量 = -Jfactor*这个点的数值*（周围六个点的数值之和）
double AvrE(int ***&a,int SIZE) {
    double E = 0;
    for (int i = 0; i < SIZE; i++)
        for (int j = 0; j < SIZE; j++)
            for (int k = 0; k < SIZE; k++)
                for (int n = 0; n < 6; n++) {
                    //计算UnitE并累加
                    int x1 = (i + dx[n] + SIZE) % SIZE, y1 = (j + dy[n] + SIZE) % SIZE, z1 = (k + dz[n] + SIZE) % SIZE;
                    E += -Jfactor * a[i][j][k] * a[x1][y1][z1];
                }
    return (double) E / (SIZE * SIZE * SIZE * 2);
}

//系统初始化 random=Ture 随机初始化 random = false 全赋相同值即全为1或全为-1
void initcond(int ***&a, bool random = true,int SIZE=32) {
    mt19937 ran(std::chrono::system_clock::now().time_since_epoch().count());
    uniform_int_distribution<int> v(0,1);
    if (random) {
        for (int i = 0; i < SIZE; i++)
            for (int j = 0; j < SIZE; j++)
                for (int k = 0; k < SIZE; k++)
                        a[i][j][k] = v(ran)*2-1;

    } else {
        int init = v(ran)*2-1;
        for (int i = 0; i < SIZE; i++)
            for (int j = 0; j < SIZE; j++)
                for (int k = 0; k < SIZE; k++)
                    a[i][j][k] = init;
    }
}

// 给定simulation模式
//mode=1 测量tem at 0~10，(mcs较小只有E,M测量值较为精确)
//mode=2 测量tem at 3~7.5 (mcs较小只有E,M测量值较为精确)
//mode=3 测量tem at 4.4~4.6 (mcs很大，U4也较为精确)
//mode=4 测量tem at 4.45~4.55 (mcs很大，通过计算U4,进一步精确相变点温度Tc，其目前精确解约为4.511523)
//温度范围可不断调试改变 与上述可能有所差异
vector<double> Tem(int mode) {
    vector<double> tem;
    if (mode == 1) {
        double Tmin = 2, Tmid0 = 3.5, Tmid1 = 4.2, Tmid2 = 4.8, Tmid3 = 6, Tmax = 8;
        int len1 = 15, len2 = 14, len3 = 18, len4 = 24, len5 = 20;
        for (int i = 0; i < len1; i++)
            tem.emplace_back((Tmid0 - Tmin) / (len1 ) * (i + 1) + Tmin);
        for (int i = 0; i < len2; i++)
            tem.emplace_back((Tmid1 - Tmid0) / len2 * i + Tmid0);
        for (int i = 0; i < len3; i++)
            tem.emplace_back((Tmid2 - Tmid1) / len3 * i + Tmid1);
        for (int i = 0; i < len4; i++)
            tem.emplace_back((Tmid3 - Tmid2) / len4 * i + Tmid2);
        for (int i = 0; i < len5; i++)
            tem.emplace_back((Tmax - Tmid3) / len5 * i + Tmid3);
    } else if (mode == 2) {
        double T1 = 2, T2 = 4.2, T3 = 4.8, T4 = 6, T5 = 7.5;
        int len1 = 16, len2 = 30, len3 = 18, len4 = 15;
        for (int i = 0; i < len1; i++)
            tem.emplace_back((T2 - T1) / len1 * i + T1);
        for (int i = 0; i < len2; i++)
            tem.emplace_back((T3 - T2) / len2 * i + T2);
        for (int i = 0; i < len3; i++)
            tem.emplace_back((T4 - T3) / len3 * i + T3);
        for (int i = 0; i < len4; i++)
            tem.emplace_back((T5 - T4) / len4 * i + T4);
    } else if (mode == 3) {
        double T1 = 4.4, T2 = 4.6;
        int len1 = 11;
        for (int i = 0; i < len1; i++)
            tem.emplace_back((T2 - T1) / (len1 - 1) * i + T1);
    } else if (mode == 4) {
        double T1 = 4.45, T2 = 4.55;
        int len1 = 11;
        for (int i = 0; i < len1; i++)
            tem.emplace_back((T2 - T1) / (len1 - 1) * i + T1);
    }else if(mode==5)
    {double T1 = 4.48, T2 = 4.54;
        int len1 = 13;
        for (int i = 0; i < len1; i++)
            tem.emplace_back((T2 - T1) / (len1 - 1) * i + T1);

    }
    else if(mode==6)
    {double T1 = 4.48, T2 = 4.51;
        int len1 = 7;
        for (int i = 0; i < len1; i++)
            tem.emplace_back((T2 - T1) / (len1 - 1) * i + T1);

    }
    else if(mode==7)
    {double T1 = 4.51, T2 = 4.54;
        int len1 = 6;
        for (int i = 0; i < len1; i++)
            tem.emplace_back((T2 - T1) / (len1) * (i+1) + T1);

    }
    else if(mode==8)
    {double T1 = 4.506, T2 = 4.516;
        int len1 = 11;
        for (int i = 0; i < len1; i++)
            tem.emplace_back((T2 - T1) / (len1 - 1) * i + T1);
    }
    else if(mode==9)
    {double T1 = 4.508, T2 = 4.514;
        int len1 = 13;
        for (int i = 0; i < len1; i++)
            tem.emplace_back((T2 - T1) / (len1-1) * i + T1);
    }

    else if(mode==10)
    {double T1 = 4.51, T2 = 4.525;
        int len1 = 3;
        for (int i = 0; i < len1; i++)
            tem.emplace_back((T2 - T1) / (len1) * i + T1);
    }
    else if(mode==11)
    {double T1 = 2, T2 = 7.5;
        int len1 = 55;
        for (int i = 0; i < len1; i++)
            tem.emplace_back((T2 - T1) / len1  * i + T1);
    }
    else if(mode==12)
    {double T1 = 8, T2 = 0;
        int len1 = 80;
        for (int i = 0; i < len1; i++)
            tem.emplace_back((T2 - T1) / len1  * i + T1);
    }


    return tem;
}

//在tembalance函数中使用，计算两个系统能量差，磁化强度差等
vector<double> var(vector<vector<vector<double>>> &v, int count)
{
    double avrdeltaE = 0, avrdeltaM = 0, e2 = 0, m2 = 0, e1 = 0, m1 = 0;
    vector<vector<vector<double>>> em;
    vector<double> res;
    em.insert(em.begin(), v.end() - count, v.end());
    for (auto i:em)
    {
        avrdeltaE += i[0][2];
        avrdeltaM += i[1][2];
        e1 += i[0][0] + i[0][1];
        m1 += i[1][0] + i[1][1];
    }
    avrdeltaE /= em.size();
    avrdeltaM /= em.size();
    e1 /= 2 * em.size();
    m1 /= 2 * em.size();
    for (auto i:em)
    {
        e2 += (i[0][2] - avrdeltaE) * (i[0][2] - avrdeltaE);
        m2 += (i[1][2] - avrdeltaM) * (i[1][2] - avrdeltaM);
    }
    e2 /= em.size() * abs(e1);
    m2 /= em.size() * abs(m1);
    res.emplace_back(e2);
    res.emplace_back(m2);
    res.emplace_back(e1);
    res.emplace_back(m1);
    res.emplace_back(abs(avrdeltaE));
    res.emplace_back(abs(avrdeltaM));
    return res;
}

//通过两个系统从不同的初始条件开始迭代，计算平衡时间，本文件主函数中未调用
vector<vector<vector<double>>> tembalance(int ***&a, int ***&b, long long int mcs, double temperature, int cond1, int cond2, int mode, int size)
{
    initcond(a, cond1, size);
    initcond(b, cond2, size);
    int count = 10;
    double e2v=1e-1,m2v=3*1e-2;
    long long int steps = 0;
    vector<double> E, M;
    vector<vector<double>> EM, avrdeltaEM;
    vector<vector<vector<double>>> RES;
    vector<double> resvar;
    if (mode == 1)
    {
        time_t t1 = clock();
        for (steps = 0; steps < 1000 * mcs; steps++)
        {
            if (steps % mcs == 0)
            {
                E.emplace_back(AvrE(a,size));
                E.emplace_back(AvrE(b,size));
                E.emplace_back(AvrE(a,size) - AvrE(b,size));
                E.emplace_back((double)steps / (size * size * size));
                M.emplace_back(AvrM(a,size));
                M.emplace_back(AvrM(b,size));
                M.emplace_back(AvrM(a,size) - AvrM(b,size));
                M.emplace_back(double(clock() - t1) / CLOCKS_PER_SEC);

                EM.emplace_back(E);
                EM.emplace_back(M);
                RES.emplace_back(EM);
                E.clear();
                M.clear();
                EM.clear();
            }

            clusterflip(a, temperature, size);
            clusterflip(b, temperature, size);
            if (steps % mcs == 0)
            {
                // cout << "clucter stpes = " << (double)steps / (size * size * size) << "  time = " << double(clock() - t1) / CLOCKS_PER_SEC << "s" << endl;
                if (steps > count * mcs)
                {
                    resvar = var(RES, count);
                    // cout << " e2= " << resvar[0] << " m2= " << resvar[1] << endl;
                    if (resvar[0] < e2v && resvar[1] < m2v)
                    {

                        resvar.emplace_back(double(clock() - t1) / CLOCKS_PER_SEC / 60);
                        cout << "Tem = " << temperature << " mode =1 "<<" steps = "<<steps
                             << " Balance time = " << double(clock() - t1) / CLOCKS_PER_SEC / 60 << "min"
                             << endl;
                        break;
                    }
                }
            }
        }
    }
    else if (mode == 2)
    {
        // mcs = mcs * size * size * size / 5;
        time_t t1 = clock();
        for (steps = 0; steps < 5000 * mcs; steps++)
        {

            if (steps % mcs == 0)
            {
                E.emplace_back(AvrE(a,size));
                E.emplace_back(AvrE(b,size));
                E.emplace_back(AvrE(a,size) - AvrE(b,size));
                E.emplace_back(steps / (size * size * size));
                M.emplace_back(AvrM(a,size));
                M.emplace_back(AvrM(b,size));
                M.emplace_back(AvrM(a,size) - AvrM(b,size));
                M.emplace_back(double(clock() - t1) / CLOCKS_PER_SEC);
                EM.emplace_back(E);
                EM.emplace_back(M);
                RES.emplace_back(EM);
                E.clear();
                M.clear();
                EM.clear();
            }

            // cout<<" a = "<<clusterflip(a, temperature, size)
            //     <<" b = "<<clusterflip(b, temperature, size)<<endl;
            clusterflip(a, temperature, size);
            clusterflip(b, temperature, size);

            if (steps % mcs == 0)
            {
                // cout << "mcsstep = " << (double)steps / (size * size * size) << "  time = " << double(clock() - t1) / CLOCKS_PER_SEC << "s " << endl;
                if (steps > count * mcs)
                {
                    resvar = var(RES, count);

                    // cout << " e2= " << resvar[0] << " m2= " << resvar[1] << endl;
                    if (resvar[4] < e2v && resvar[5] < m2v)
                    {
                        resvar.emplace_back(steps);
                        resvar.emplace_back(double(clock() - t1) / CLOCKS_PER_SEC);
                        cout << "sys " << size << " under temperture at " << temperature << " at  mode = " << mode
                             << " balance time is " << double(clock() - t1) / CLOCKS_PER_SEC << "s"
                             << " steps =  " << steps <<" deltaE = "<<resvar[4]<<" deltaM = "<<resvar[5]
                             << endl;
                        break;
                    }
                }
            }
        }
    }
    else if (mode == 3)
    {
        // mcs = mcs * size * size * size / 5;
        time_t t1 = clock();
        int i=0;
        for (steps = 0; steps < mcs*350; steps++)
        {

            if (steps % mcs == 0)
            {
                E.emplace_back(AvrE(a,size));
                E.emplace_back(AvrE(b,size));
                E.emplace_back(AvrE(a,size) - AvrE(b,size));
                E.emplace_back(steps / (size * size * size));
                M.emplace_back(AvrM(a,size));
                M.emplace_back(AvrM(b,size));
                M.emplace_back(AvrM(a,size) - AvrM(b,size));
                M.emplace_back(double(clock() - t1) / CLOCKS_PER_SEC);
                EM.emplace_back(E);
                EM.emplace_back(M);
                RES.emplace_back(EM);
                E.clear();
                M.clear();
                EM.clear();
                RES[i++][0].emplace_back(steps);
            }

            clusterflip(a, temperature, size);
            clusterflip(b, temperature, size);
        }
    }
    if(mode!=3)
        RES[0].emplace_back(resvar);
    return RES;
}

//对于给定温度 给定大小为Size系统达到平衡的步长mcs 给定系统达到平衡后的采样间隔interval 给定初始化条件随机赋值or不随机赋值
//返回采样结束后获得的系统的一些参数如平均能量AvrE，平均磁化强度AvrM，比热容C,磁化率X,磁化累积量U4,通过精确计算U4可得到相变点温度Tc
vector<double> start(int ***&a, int SIZE,double temperature, long long int mcs, int interval, bool init = true) {
    //初始化特定温度下的系统 initcond=Ture 随机初始化 initcond = false 全赋相同值即全为1或全为-1
    initcond(a, init,SIZE);
//    double pdd; //接受率
    pdd = 1 - exp(-2 * Jfactor / temperature);
//    mcs/=(pdd*pdd);
//    interval/=(pdd*pdd);
    if(temperature>4.5)
        mcs*=5;
    vector<vector<double>> res(6);
    vector<double> RES;
    for (int step = 0; step < 100 * mcs; step++) {
        // time_t t1 = clock();
        clusterflip(a, temperature,SIZE);

        if (step >= 30 * mcs && step % interval == 0) {
            double ae = AvrE(a,SIZE);
            double am = AvrM(a,SIZE);
            res[0].emplace_back(ae);
            res[1].emplace_back(am);
            res[2].emplace_back(am * am);
            res[3].emplace_back(am * am * am * am);
            res[4].emplace_back(ae * ae);
            res[5].emplace_back(ae * ae * ae * ae);
        }
        // cout<<"T0 "<<(double)(clock()-t1)/CLOCKS_PER_SEC<<"s."<<endl;
    }
    int LEN = res[0].size();
    double AE = 0, AM = 0, M2 = 0, M4 = 0, E2 = 0, E4 = 0, Evar = 0;

    for (int i = 0; i < LEN; i++) {
        AE += res[0][i] / LEN;
        AM += res[1][i] / LEN;
        M2 += res[2][i] / LEN;
        M4 += res[3][i] / LEN;
        E2 += res[4][i] / LEN;
        E4 += res[5][i] / LEN;
    }
    for (int i = 0; i < LEN; i++)
        Evar += (res[0][i] - AE) * (res[0][i] - AE) / LEN;

    RES.emplace_back(AE);                                                                //平均能量
    RES.emplace_back(AM);                                                                //平均磁化强度
    RES.emplace_back(1 - M4 / (3 * M2 * M2));                                            //U4阶矩
    RES.emplace_back(SIZE*SIZE*SIZE * Evar / (temperature * temperature));           //比热 = 能量的方差/温度平方
    RES.emplace_back(SIZE*SIZE*SIZE * (E2 - AE * AE) / (temperature * temperature)); //比热' 方差公式<(E-<E>)^2>=<E^2>-<E>^2
    RES.emplace_back((M2 - AM * AM) / temperature);                                      //磁化率 <(M-<M>)^2>/tem=(<M^2>-<M>^2)/tem

    return RES;
}

void init_a(int ***&a, int SIZE) {
    a = new int **[SIZE];
    for (int i = 0; i < SIZE; i++) {
        a[i] = new int *[SIZE];
        for (int j = 0; j < SIZE; j++) {
            a[i][j] = new int[SIZE];
        }
    }
}

int main(int argc, char **argv) {

    //Tc=4.511523
    char *ptr;
    int SIZE = strtol(argv[1],&ptr,10);
    init_a(canvas, SIZE);                                   //分配空间
    long long int mcs = strtol(argv[2],&ptr,10);    //假设达到稳态所用的step,每个温度一共循环100mcs
    long int interval = strtol(argv[3],&ptr,10);  //达到稳态后（为了更准确，过了20个mcs后进行采样)采样step间隔S
    bool init = false;                                         //系统初始化条件
    int mode = strtol(argv[4],&ptr,10);          //测试特定模式  mode1,2仅仅为了得到E,M mcs可小一点，mode 3,4为了得到精确的U4 mcs需较大
    int spec = strtol(argv[5],&ptr,10);             //同一参数下的第几个文件
    ofstream file1;
    ostringstream path(R"(D:\res_)", ios_base::ate);
    path << SIZE << "_" << mcs << "_" << interval << "_" << mode <<"_"<<spec<< ".csv";
    file1.open(path.str(), ios::out | ios::trunc);
    cout << "size= " << SIZE << " mcs= " << mcs << " interval =" << interval << " mode= " << mode << endl;

    vector<vector<double>> RES;
    vector<double> tem;
    tem = Tem(mode);
    int tot_len = tem.size();
    for (int i = 0; i < tot_len; i++) {
        time_t t1 = clock();
        RES.emplace_back(start(canvas,SIZE, tem[i], mcs*(1+(tem[i]/Tc)), interval*(1+(tem[i]/Tc)), init));
        cout << "Temperature = " << tem[i] << "   AverageEnergy = " << RES[i][0] << "   AverageMag = " << RES[i][1]
             << "   U4 = " << RES[i][2] << "   Cv1 = " << RES[i][3] << "   Cv2 = " << RES[i][4] << "   X = "
             << RES[i][5];
        time_t t2 = clock();
        RES[i].emplace_back((double) (t2 - t1) / CLOCKS_PER_SEC);
        cout << "   single temperature costs " << RES[i][6] / 60 << "min." << endl;

        //输出文件 第一列为温度，第二列为能量，第三列为磁化强度，第四列为U4，第五列第六列为两种方法计算的比热，第六列为磁化率
        file1 << tem[i] << ",";

        for (int j = 0; j < RES[0].size(); j++)

            file1 << RES[i][j] << ",";
        file1 << endl;
    }
    file1.close();

    double tot_time;
    for (int i = 0; i < tot_len; i++)
        tot_time += RES[i][6];

    cout << " Complete!" << endl;
    cout << " Total cost " << tot_time / 60 << "min" << endl;

//    getchar();
    return 0;
}
