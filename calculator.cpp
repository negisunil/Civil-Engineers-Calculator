#include<bits/stdc++.h>
#include<iostream>
#include<math.h>
#include<string>
#include<map>
#include<vector>
#include<conio.h>

#define EI 90000

using namespace:: std;


////////////////////////////////////////////////////////////SOM Functions//////////////////////////////////////////////////////////////////////////

void SlopeDef_SSB(float span, float P)
{
    char ch;
    cout<<" Please choose the type of loading on steel simply supported beam :"<<endl;
    cout<<" press a->for central point load\n  b->for equal and unlike point moment at ends"<<endl;
    cout<<" press c->for central point moment\n press d->for UDL\n press e->for UVL with may at centre\n f->UVL with may at support";
    cin>>ch;
    float slope_s,slope_m,def_m;
    slope_s=slope_m=def_m=0;
    
    switch(ch)
    {
        case 'a': slope_s=P*pow(span,2)/(16*EI);
                  def_m= P*pow(span,3)*1000/(48*EI);
                  break;
                  
        case 'b': slope_s= P*span/(24*EI);
                  slope_m= slope_s*2;
                  break;
        
        case 'c':  slope_s=P*span/(2*EI);
                   def_m=P*pow(span,2)*1000/(8*EI);
                   break;
                   
        case 'd': slope_s=P*pow(span,3)/(24*EI);
                  def_m= P*pow(span,4)*5000/(384*EI);
                  break;   
                  
        case 'e': slope_s=P*pow(span,3)*5/(192*EI);
                  def_m= P*pow(span,4)*1000/(120*EI);
                  break;
        
        case 'f': slope_s=P*pow(span,3)*5/(192*EI);
                  def_m= P*pow(span,4)*1000/(120*EI);
                  break;
                  
        default : cout<<" invalid input";
                   break;
                  
    }
    
    cout<<"\n The slope at supports and mid section of beam is :"<<slope_s<<" and "<<slope_m<<" (in radians)";
    cout<<"\n deflection at mid section is "<<def_m<<" mm";
    return;
}

void SlopeDef_Cantilever(float span, float P)
{
    char ch;
    cout<<" Please choose the type of loading on steel cantilever beam :"<<endl;
    cout<<" press a->for point load at free end\n  b->for point moment at free end"<<endl;
    cout<<" press c->for UDL\n press d->for UVL with max at support";
    cin>>ch;
    float slope_s,slope_e,def_e;
    slope_s=slope_e=def_e=0;
    
    switch(ch)
    {
        case 'a': slope_e=P*pow(span,2)/(2*EI);
                  def_e= P*pow(span,3)*1000/(3*EI);
                  break;
                  
        case 'b': slope_e= P*span/EI;
                  def_e=P*pow(span,2)/(2*EI);
                  break;
        
        case 'c':  slope_e=P*pow(span,3)/(6*EI);
                   def_e=P*pow(span,4)*1000/(8*EI);
                   break;
                   
        case 'd': slope_e=P*pow(span,3)/(24*EI);
                  def_e= P*pow(span,4)*1000/(30*EI);
                  break;   
                  
        default : cout<<" invalid input";
                   break;
                  
    }
    
    cout<<"\n The slope at support and free end of beam is :"<<slope_s<<" and "<<slope_e<<" (in radians)";
    cout<<"\n deflection at free end is "<<def_e<<" mm";
    return;
}

////////////////////////////////////////////////////////////Structure Functions//////////////////////////////////////////////////////////////////////////

void FEM(float span, float W)
{
    char ch;
    cout<<" Please choose the type of loading on steel Fixed beam :"<<endl;
    cout<<" press a->for point load at x distance from left support "<<endl;
    cout<<" press b->for point moment at x distance from left support"<<endl;
    cout<<" press c->for UDL entire span\n press d->for UVL with max at left support"<<endl;
    cout<<" press e->for UVL with max at right support\n press f->for UVL with max at center"<<endl;
    cin>>ch;
    float M1,M2,x;
   
    switch(ch)
    {
        case 'a': cout<<" \n Enter the value of x in meters";
                  cin>>x;
                  M1=W*x*pow(span-x,2)/pow(span,2);
                  M2=W*(span-x)*pow(x,2)/pow(span,2);
                  cout<<"\n left and right supports will have anticlockwise and clockwise moment respectively "<<endl;
                  break;
                  
        case 'b': cout<<" \n Enter the value of x in meters";
                  cin>>x;
                  M1=W*(span-x)*(3*x-span)/pow(span,2);
                  M2=W*x*(2*span-3*x)/pow(span,2);
                  cout<<"\n left and right supports will have moment of same nature as applied"<<endl;
                  break;
        
        case 'c': M1=M2=W*pow(span,2)/12;
                  cout<<"\n left and right supports will have anticlockwise and clockwise moment respectively "<<endl;
                  break;
                   
        case 'd': M1=W*pow(span,2)/20;
                  M2=W*pow(span,2)/30;
                  cout<<"\n left and right supports will have anticlockwise and clockwise moment respectively "<<endl;
                  break;  
                  
        case 'e': M1=W*pow(span,2)/30;
                  M2=W*pow(span,2)/20;
                  cout<<"\n left and right supports will have anticlockwise and clockwise moment respectively "<<endl;
                  break;  
                  
        case 'f': M1=M2=W*pow(span,2)*5/96;
                  cout<<"\n left and right supports will have anticlockwise and clockwise moment respectively "<<endl;
                  break;           
                  
        default : cout<<" invalid input";
                   break;
                  
    }
    
    cout<<"\n Moments at left and right supports are :"<<M1<<" and "<<M2<<" (in KN-m) repectively";
    return;
}

///////////////////////////////////////////////////////////Environment Functions//////////////////////////////////////////////////////////////////////

void BOD(float t, float k20, float T , float BODu)
{
    float k=k20*pow(1.047,T-20);
    float BODt=BODu*(1-exp(-k*t));
    float BODr=BODu*exp(-k*t);
    cout<<"\n BOD consumed and remaining after "<<t<<" days at "<<T<<" degree Celcius are "<<BODt<<" and "<<BODr;
}

void settlingvel(float d) //input dia is in mm
{
    double Vs,nu=0.000001; //assuming inorganic solids G=2.65 and temperature =20 degree C 
    
    if(d<0.1)
    {
        cout<<"\n Laminar flow";
        Vs=9.81*1.65*pow(d/1000,2)/(18*nu);
    }
    
    else if(d>0.1 && d<1)
    {
        cout<<"\n Transition flow";
        Vs=pow(9.81*1.65*pow(d/1000,1.6)/(13.88*pow(nu,0.6)),0.714);
    }
    
    else
    {
        cout<<"\n Turbulent flow";
        Vs=1.8*sqrt(9.81*d*1.65/1000);  //Newton's equation
    }
    
    cout<<"\n Settling velocity of the inorganic solid is "<<Vs<<" m/s";
    return;
}

void hardness_alk(float Ca2, float Mg2, float HCO3, float CO3) //all input in mg/l
{
    float TH, CH, NCH, alk;
    TH=((Ca2/20)+(Mg2/12))*50 ;
    alk=((HCO3/61)+(CO3/30))*50;
    
    if(TH>alk)
    CH=alk;
    
    else
    CH=TH;
    
    NCH=TH-CH;
    
    cout<<"\n Total, carbonate & non-carbonate hardness of given sample are "<<TH<<", "<<CH<<" and "<<NCH<<"(all in mg/L)";
    cout<<"\n Alkalinity of the sample is "<<alk<<" mg/L";
    
    if(TH<=55)
    cout<<"\n Soft water ";
    
    else if(TH>=56 &&TH<=100)
    cout<<"\n Slightly hard water";
    
    else if(TH>=101 &&TH<=200)
    cout<<"\n Moderately hard water";
    
    else
    cout<<"\n Very hard water";
    
    return;
}

////////////////////////////////////////////////////////////Numerical Methods//////////////////////////////////////////////////////////////////////

//NEWTON RAPHSON
double fx(double x)
{
    return 3*x +sin(x)-pow(2.71828,x);
}
 
double derivative(double x)
{
    return 3 +cos(x)-pow(2.71828,x);
}
 
void newtonRaphson(double x, double acc)
{
    double h = fx(x) / derivative(x);
    while (abs(h) >= acc)
    {
        h = fx(x)/derivative(x);
  
        // x(i+1) = x(i) - f(x) / f'(x) 
        x = x - h;
    }
 
    cout << "The value of the root is : " << x;
}

//Bisection method for cubic equation

double p,q,r,s;   //Global Variables

double f(double x)
{
    return p*x*x*x +q*x*x +r*x + s;
}
 

void bisection(double a, double b, double acc)
{
    if (f(a) * f(b) >= 0)
    {
        cout << "You have not assumed right a and b\n";
        return;
    }
 
    double c = a;
    while ((b-a) >= acc)
    {
    
        c = (a+b)/2;

        if (f(c) == 0.0)
            break;
 
        else if (f(c)*f(a) < 0)
            b = c;
        else
            a = c;
    }
    cout << "The value of root is : " << c;
}

///////////////////////////////////////////////////////////Statistics functions///////////////////////////////////////////////////////////////////

float mean_ungrouped(  vector<float> v)
{
	float cnt,sum;
	cnt=sum=0;
	for(int i=0; i<v.size();i++)
	{
		sum+=v[i];
		cnt++;
	}
	return (sum/cnt);
}


float variance(vector<float> v)
{
	float avg=mean_ungrouped(v);
	float sum,cnt;
	sum=cnt=0;
	for(int i=0; i<v.size(); i++)
	{
		sum+=(avg-v[i])*(avg-v[i]);
		cnt++;
	}
	
	return (sum/(cnt-1));	
}

float std_dev(vector<float> v)
{
	float res=sqrt(variance(v));
	return res;
}

float coeff_variance(vector<float> v)
{
	float sd=std_dev(v);
	float avg=mean_ungrouped(v);
	
	return (sd/avg);
}

float mode_ungrouped(vector<float> v)
{
	map<float,int> ump;
	
	for(int i=0; i<v.size(); i++ )
	ump[v[i]]++;
	map<float,int> ::iterator it;
	int m=0;
	float mo;
	
	for(it=ump.begin(); it!=ump.end(); it++)
	{
		if(it->second>m)
		{ 
		    m=it->second;
		    mo=it->first;
		}
		
	}
	
	return mo;
}

float median_ungrouped(vector<float> v)
{
	sort(v.begin(),v.end());
	int n=v.size();
	float median;
	if(n%2!=0)
	median=(v[n/2] +v[(n/2) -1])/2 ;
	
	else
	median=v[n/2];
	
	return median;
	
}

////////////////////////////////////////////////////////////////Basic Mathematical Functions///////////////////////////////////////////////////////

int factorial(int n)
{
	
	if(n==0||n==1)
	return 1;
	
	else
	return n*factorial(n-1);
}

float combination(int n, int r)
{
	float nCr=factorial(n)/(factorial(n-r) *factorial(r));
	
	return nCr;
	
}

float permutaion(int n, int r)
{
	float nPr=factorial(n)/factorial(n-r);
	
	return nPr;
}
void swap(double &x, double &y)
{
	x=x+y;
	y=x-y;
	x=x-y;
}

/////////////////////////////////////////////////////////////////////Equations////////////////////////////////////////////////////////////////////

void roots_quadratic(double a, double b, double c)
{
    if(a==0)
    {
        cout<<"\t invalid input";
        return;
    }
    
    double D=b*b-4*a*c ;
    double rootD=sqrt(abs(D));
    
    if(D==0)
    {
        cout<<"\troots are real and equal"<<endl;
        cout<<"\tvalue of both the roots is: "<<(-b/(2*a));
        return;
    }
    
    if(D>0)
    {
        double root1,root2;
        cout<<"\troots are real and distinct"<<endl;
        root1=(-b+rootD)/(2*a);
        root2=(-b-rootD)/(2*a);
        cout<<"\tvalues of the roots are: "<<root1<<"\t"<<root2;
        return;
    }
    if(D<0)
    {
        cout<<"\troots are imaginary"<<endl;
        double re,im;
        re=-b/(2*a);
        im= rootD/(2*a);
        cout<<"\tvalues of the roots are: "<<re<<" + i"<<im<<"\t "<<re<<" - i"<<im;
        return;
    }
    
    return;
}

void Lequation_2var(double coeff[2][3])
{
    if((coeff[0][0]/coeff[1][0])==(coeff[0][1]/coeff[1][1]) && (coeff[0][1]/coeff[1][1])!=(coeff[0][2]/coeff[1][2]))
    {
        cout<<"\t No Solution Possible of this system of equations"<<endl;
        return;
    }
    else if((coeff[0][0]/coeff[1][0])==(coeff[0][1]/coeff[1][1]) && (coeff[0][1]/coeff[1][1])==(coeff[0][2]/coeff[1][2]))
    {
        cout<<"\t Infinitely many Solutions Possible of this system of equation"<<endl;
        return;
    }
    else
    {
        //CRAMER'S Rule
        double D,d1,d2;
        D=coeff[0][0]*coeff[1][1] - coeff[1][0]*coeff[0][1];
        d1=coeff[0][2]*coeff[1][1] - coeff[1][2]*coeff[0][1];
        d2=coeff[0][0]*coeff[1][2] - coeff[1][0]*coeff[0][2];
        
        cout<<"\t Unique Solution of given system is: y="<<d1/D<<"  y="<<d2/D;
        return;
    }
    
}
//declare det_3x3 before it
double det_3x3(double mat[3][3]);

void Lequation_3var(double coeff[3][4])
{
    double D,d1,d2,d3,i,j;
    double d[3][3]={ {coeff[0][0],coeff[0][1],coeff[0][2]},
                  {coeff[1][0],coeff[1][1],coeff[1][2]},
                  {coeff[2][0],coeff[2][1],coeff[2][2]} };
                  
    double a[3][3]={ {coeff[0][3],coeff[0][1],coeff[0][2]},
                  {coeff[1][3],coeff[1][1],coeff[1][2]},
                  {coeff[2][3],coeff[2][1],coeff[2][2]} };  
                  
    double b[3][3]={ {coeff[0][0],coeff[0][3],coeff[0][2]},
                  {coeff[1][0],coeff[1][3],coeff[1][2]},
                  {coeff[2][0],coeff[2][3],coeff[2][2]} };   
                  
    double c[3][3]={ {coeff[0][0],coeff[0][1],coeff[0][3]},
                  {coeff[1][0],coeff[1][1],coeff[1][3]},
                  {coeff[2][0],coeff[2][1],coeff[2][3]} };  
                  
                  
    D = det_3x3(d);
    d1= det_3x3(a);
    d2= det_3x3(b);
    d3= det_3x3(c);
    
    if(!D)
    {
        if(d1==0 || d2==0 || d3==0)
        {
            cout<<"\t Infinitely many Solutions Possible"<<endl;
            return;
        }
        
        else
        {
            cout<<"\t No Solution Possible of this system"<<endl;
            return;
        }
    }
    
    else
    {
        cout<<"\t Unique Solution of given system is: x="<<d1/D<<"  y="<<d2/D<<"  z="<<d3/D;
        return;
    }
    
}


////////////////////////////////////////////////////////////////MATRICES////////////////////////////////////////////////////////////////////

double det_3x3(double mat[3][3])
{
    double ans;
    ans=mat[0][0]*(mat[1][1]*mat[2][2]-mat[1][2]*mat[2][1]) 
        - mat[0][1]*(mat[1][0]*mat[2][2]-mat[2][0]*mat[1][2]) 
        + mat[0][2]*(mat[1][0]*mat[2][1]-mat[2][0]*mat[1][1]);
    return ans;     
}

double det_2x2(double mat[2][2])
{
    double ans;
    ans=mat[0][0]*mat[1][1]-mat[1][0]*mat[0][1];
    return ans;     
}
void eigenvalues_2x2(double  mat[2][2])
{
    double a,b,c;
    a=1;
    b=-(mat[0][0]+mat[1][1]);
    c=mat[0][0] *mat[1][1] - mat[0][1] *mat[1][0] ;
    
    roots_quadratic(a,b,c);
    
}

void adjoint_2x2(double mat[2][2])
{
	swap(mat[0][0],mat[1][1]);
	mat[0][1]*=-1;
	mat[1][0]*=-1;
	cout<<"\n adjoint is : \n";
	cout<<mat[0][0]<<"  "<<mat[0][1]<<endl;
	cout<<mat[1][0]<<"  "<<mat[1][1];
}

void transpose_3x3(double mat[3][3])
{
    int i,j;
    for(i=0; i<3 ;i++)
    {
        for(j=i+1;j<3;j++)
        {
            mat[i][j] =mat[i][j]+mat[j][i];
            mat[j][i]=mat[i][j]-mat[j][i];
            mat[i][j]=mat[i][j]-mat[j][i];
        }
    }
    
    for(i=0; i<3 ;i++)
    {
        for(j=0;j<3;j++)
        cout<<"\t"<<mat[i][j];
        
        cout<<endl;
    }

}

///////////////////////////////////////////////////////////////////Probability/////////////////////////////////////////////////

void bionomial_dist(float p,int n, int r)
{
	if(p>1||p<0||r>n||r<0||n<0)
	{
		cout<<" wrong inputs";
		return;
	}
	float probability=combination(n,r) *pow(p,r)*pow(1-p,n-r);
	
	cout<<"\t Probability of the given event is: "<<probability;
}

void poisson_dist(float x, float mean)
{
    if(mean<0||r<0)
	{
		cout<<" wrong inputs";
		return;
	}
	float probability;
	
	probability=pow(mean,x)*exp(-mean)/factorial(x);
	
	cout<<"\t Probability of the given event is: "<<probability;
}


int main()
{
	system("color 0C");
	cout<<"\n\n#################################################### WELCOME TO CIVIL ENGINEERS CALCULATOR #################################################\n\n";
	int in;
	char ch;
	cout<<" Choose the function to be performed : "<<endl;
	cout<<" press 1 for SOM functions\n press 2 for structre analysis function"<<endl;
	cout<<" press 3 for environmental engineering functions\n press 4 for numerical methods\n";
	cout<<" press 5 for equations\n press 6 for probability distributions"<<endl;
	cout<<" press 7 for matrices functions \n press 8 for statistical functions";
	cin>>in;
	
	switch(in)
	{
		case 1: cout<<" Enter the beam type either SSB or cantilever, length in meters and load or moment in KN or KN-m\n";
		        char beamtype[20];
		        float l,P;
		        cin>>beamtype>>l>>P;
		        if(beamtype=="SSB")
		        SlopeDef_SSB(l,P);
		        else if(beamtype=="cantilever")
		        SlopeDef_Cantilever(l,P);
		        else
		        cout<<"wrong input";
		        break;
		        
		case 2: cout<<" Enter the fixed beam length in meters and load or moment in KN or KN-m\n";
		        float span,W;
		        cin>>span>>W;
		        FEM(span,W);
		        break;
				
		case 3: char ch;
		        cout<<" select the function to be performed:\n press a for BOD\n";
		        cout<<" press b for settling velocity\n press c for alkalinity and hardness\n";
				cin>>ch;
				if(ch=='a')
				{
					float t,T,Bu,k;
					cout<<"Enter K(base e at 20 degree C) per day, current temperature in degree C,";
					cout<<" ultimate BOD in mg/L and time after which bod required in days";
					cin>>k>>T>>Bu>>t;
					BOD(t,k,T,Bu);	
						   }
				else if(ch=='b')
				{
					float d;
					cout<<"Enter diameter of particles in mm\n";
					cin>>d;
					settlingvel(d);
				}
				else if(ch=='c')
				{
					float Ca,Mg,HCO3,CO3;
					cout<<"Enter the concentration of Ca, Mg, HCO3- and CO3-- in mg/L\n";
					cin>>Ca>>Mg>>HCO3>>CO3;
					hardness_alk(Ca,Mg,HCO3,CO3);
				}
				else
				cout<<"wrong input";
				break;	
				
		case 4: char nm;
		        cout<<" select the method to be performed:\n press a for Newton Raphson \n";
		        cout<<" press b for Bisection\n";
				cin>>nm;
				double accuracy;			   		       
		        if(nm=='a')
				{
					cout<<"This method is only availabe for 3x+sin(x)-e^x type equation\n";
					double x0;
                    cout<<"Enter the initial guess and required accuracy\n";
                    cin>>x0>>accuracy;
                    newtonRaphson(x0,accuracy);
						   }
				else if(nm=='b')
				{
				    cout<<"This method is availabe for ax^3 + bx^2 + cx + d type cubic equation\n";
					double a,b;
                    cout<<"Enter the coefficients of cubic equation\n";
                    cin>>p>>q>>r>>s;
                    cout<<"enter the initial values and accuracy required";
                    cin>>a>>b>>accuracy;
                    bisection(a, b,accuracy);	
				}
				else
				cout<<"wrong input";
				break;
				
		case 5: char eq;
		        cout<<" select the equation type:\n press a for quadratic equation ax^2 + bx + c=0 \n";
		        cout<<" press b for linear equations in 2 variables a1x + b1y=c1 and a2x + b2y =c2\n";
		        cout<<" press c for linear equations in 3 variables\n";
				cin>>eq;
				
				if(eq=='a')
				{
				 double a,b,c;
				 cout<<"Enter the coefficients\n";
				 cin>>a>>b>>c;
				 roots_quadratic(a,b,c);
						}
				
				else if(eq=='b')
				{
					double coeffi[2][3];
					cout<<"Enter the coefficients in this order a1,b1,c1,a2,b2 and c2\n";
					for(int i=0; i<2; i++)
					{
						for(int j=0; j<3 ;j++)
						  cin>>coeffi[i][j];
					 }
					 Lequation_2var(coeffi); 
				}
				else if(eq=='c')
				{
					double coef[3][4];
					cout<<"Enter the coefficients in this order a1,b1,c1,d1,a2,b2,c2,d2,a3,b3,c3 and d3\n";
					for(int i=0; i<3; i++)
					{
						for(int j=0; j<4 ;j++)
						  cin>>coef[i][j];
					 }
					 Lequation_3var(coef); 
				}
				else
				cout<<"wrong input";
				break;
		
		case 6: char pd;
		        cout<<" select the probability distribution to be used:\n press a for Bionomial Distribution \n";
		        cout<<" press b for Poisson Distribution\n";
				cin>>pd;
				float p;								
		        if(pd=='a')
		        {
		        	int n,r;
		        	cout<<"Enter the number of trials, probability of success & no. of success for which probability is required\n";
		        	cin>>n>>p>>r;
		        	bionomial_dist(p,n,r);
				}
				else if(pd=='b')
		        {
		        	int r;
		        	cout<<"Enter the mean of distribution and no. of success for which probability is required\n";
		        	cin>>p>>r;
		        	poisson_dist(r,p);
				}
				else
				cout<<"wrong inputs";
				break;
				
		case 7: char ma;
		        cout<<" select the matrix function to be used:\n press a for 2x2 matrix determinant \n";
		        cout<<" press b for 3x3 determinant\n press c for 2x2 adjoint\n press d for 3x3 transpose\n";
		        cout<<" press e for 2x2 eigen values\n";
				cin>>ma;
				double matrix[2][2], mat[3][3];
				int i,j;
				if(ma=='a')
				{
					double d2;
					cout<<"Enter the matrix\n";
					for(i=0; i<2 ;i++)
					{
						for(j=0; j<2 ;j++)
						cin>>matrix[i][j];
					}
					d2=det_2x2(matrix);
					cout<<"\n determinant is "<<d2;
								}				
				else if(ma=='b')
				{
					double d3;
					cout<<"Enter the matrix\n";
					for(i=0; i<3 ;i++)
					{
						for(j=0; j<3 ;j++)
						cin>>mat[i][j];
					}
					d3=det_3x3(mat);
					cout<<"\n determinant is "<<d3;
								}
				else if(ma=='c')
				{
					cout<<"Enter the matrix\n";
					for(i=0; i<2 ;i++)
					{
						for(j=0; j<2 ;j++)
						cin>>matrix[i][j];
					}
					
					adjoint_2x2(matrix);
								}
			    else if(ma=='d')
				{
					cout<<"Enter the matrix\n";
					for(i=0; i<3 ;i++)
					{
						for(j=0; j<3 ;j++)
						cin>>mat[i][j];
					}
					transpose_3x3(mat);
								}
				else if(ma=='e')
				{
					cout<<"Enter the matrix\n";
					for(i=0; i<2 ;i++)
					{
						for(j=0; j<2 ;j++)
						cin>>matrix[i][j];
					}
					
					eigenvalues_2x2(matrix);
				}
				else
				cout<<"wrong inputs";
				break;
				
		case 8: vector<float> v;
		        int N;
		        float mean,mode,median,var,covar,stddev;
		        cout<<"Enter the no. of observations\n";
				cin>>N;
				float data[N];
				cout<<"\nEnter the observations now:\n ";
				
				for(int i=0; i<N;i++)
				cin>>data[i];
				
				for(int i=0; i<N;i++)
				v.push_back(data[i]);
				mean=mean_ungrouped(v);
				mode=mode_ungrouped(v);
				median=median_ungrouped(v);
				var=variance(v);
				stddev=std_dev(v);
				covar=coeff_variance(v);
				cout<<"mean = "<<mean<<" mode = "<<mode<<" median = "<<median;
				cout<<" variance = "<<var<<" std_deviation = "<<stddev<<" cv ="<<covar;
				break;
					
	}
 
    return 0;
}

