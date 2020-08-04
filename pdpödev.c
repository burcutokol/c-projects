#include <stdio.h> 
#include <stdlib.h> 
#include <string.h> 
#include <math.h> 
 
typedef struct{ 
 double x; 
 double y; 
 double z; 
} Vector; 
 
void print_vector(const Vector v)
{
	printf("(%.2f, %.2f, %.2f) \n",v.x,v.y,v.z);
	
} 
Vector sum(const Vector v1, const Vector v2)
{
	Vector summary={(v1.x+v2.x),(v1.y+v2.y),(v1.z+v2.z)};
	return summary;
}
	
Vector diff(const Vector v1, const Vector v2)
{
	Vector difference={ (v1.x-v2.x),(v1.y-v2.y),(v1.z-v2.z)};
	return difference;
} 
double dot_product(const Vector v1, const Vector v2)
{
	double dotP=(((v1.x)*(  v2.x))+((v1.y)*(v2.y))+((v1.z)*(v2.z)));
	return dotP;
}
Vector cross_product(const Vector v1, const Vector v2)
{
	Vector crossP={ ((v1.y*v2.z)-(v1.z*v2.y)),((v1.z*v2.x)-(v1.x*v2.z)),((v1.x*v2.y)-(v1.y*v2.x)) };
	return crossP;
}
double norm(const Vector v)
{
	double _norm=sqrt(pow(v.x,2)+pow(v.y,2)+pow(v.z,2));
	return _norm;
}
int is_unitvector(const Vector v)
{
	int size=norm(v);
	int unitV=0;
	if(size==1)
	{
		unitV=1;
	}
	else
	{
		unitV=0;
	}
	return unitV;
}
		
	
Vector unit(const Vector v)
{
	double size=norm(v);
	Vector unitV={0,0,0};
	if(size==0)
	{
	    printf("bos vektor...\n");
    }
    else if(size>0)
    {
        unitV.x=(v.x/size);
        unitV.y=(v.y/size);
        unitV.z=(v.z/size);
    }
    return unitV;
}
Vector multiplyby_scalar(const Vector v1, const double c)
{
	Vector newVector={(c*(v1.x)),(c*(v1.y)),(c*(v1.z))};
	return newVector;
}
double angle(const Vector v1, const Vector v2)
{
	double sizeOfV1=norm(v1);
	double sizeOfV2=norm(v2);
	double _angle=0.00;
	double dotp=dot_product(v1,v2);
	if(sizeOfV1!=0 && sizeOfV2!=0)
	{
		_angle=acos(dotp/(sizeOfV1*sizeOfV2));
	}
	else
	printf("vektorlerden biri bos vektor...\n");
	
	return _angle;	
} 
double distance(const Vector v1, const Vector v2)
{
	double _distance=sqrt(pow(v1.x-v2.x,2)+pow(v1.y-v2.y,2)+pow(v1.z-v2.z,2));
	return _distance;
}
	
int are_linearly_independent(const Vector v1, const Vector v2, const Vector v3)
{
	double sum1=((v1.x*v2.y*v3.z)+(v1.y*v2.z*v3.x)+(v1.z*v2.x*v3.y));
	double sum2=((v3.x*v2.y*v1.z)+(v3.y*v2.z*v1.x)+(v3.z*v2.x*v1.y));
	int areL=0;
	if((sum1-sum2)!=0)
	{
		areL=1;
	}
	else
	{
		areL=0;
	}
	return areL;
}
	
int are_orthogonal(const Vector v1, const Vector v2, const Vector v3)
{
	int areO=0;
	double dotp13=dot_product(v1,v3);
	double dotp23=dot_product(v2,v3);
	double dotp12=dot_product(v1,v2); 
	if(dotp13==0 && dotp23==0 && dotp12==0)
	{
		areO=1;
	}
	else
	{
		areO=0;
	}
	return areO;
}
int are_orthonormal(const Vector v1, const Vector v2, const Vector v3)
{
	int areON=0;
	int areOG=are_orthogonal(v1,v2,v3);
	int unit1=is_unitvector(v1);
	int unit2=is_unitvector(v2);
	int unit3=is_unitvector(v3);
	if(areOG==1 && unit1==1 && unit2==1 && unit3==1)
	{
		areON=1;
	}
	else
	areON=0;
	
	return areON;
}
Vector projection(const Vector v1, const Vector v2)
{
	double dotP=dot_product(v1,v2);
	double size=norm(v2);
	Vector newV={0, 0, 0};
	if(size!=0)
	{
	    double k=dotP/pow(size,2);
	    newV.x=(v2.x*k);
	    newV.y=(v2.y*k);
	    newV.z=(v2.z*k);
	}
	else
	printf("vektorlerden biri bos vektor... \nsifira bolunme hatasi...\n");
	
	return newV;
}
	 
Vector orthogonal_projection(const Vector v1, const Vector v2)
{
	Vector newV1=projection(v1,v2);
	Vector newV2={0,0,0};
	if(norm(newV1)!=0)
	{
	    newV2.x=(v1.x)-(newV1.x);
	    newV2.y=(v1.y)-(newV1.y);
	    newV2.z=(v1.z)-(newV1.z);
	}
	
	return newV2;
}
	 
int convert_2_orthogonal_basis(Vector *v1, Vector *v2, Vector *v3)
{
	int result=are_linearly_independent(*v1,*v2,*v3);
	int a=0;
	if(result==1)
	{
		Vector newV1={(*v1).x,(*v1).y,(*v1).z};
		if(dot_product(newV1,newV1)!=0)
		{
			double res=(dot_product(*v2,newV1))/(dot_product(newV1,newV1));
			Vector newV2=diff(*v2,multiplyby_scalar(newV1,res));
			double res1=(dot_product(*v3,newV1))/(dot_product(newV1,newV1));
			double res2=(dot_product(*v3,newV2))/(dot_product(newV2,newV2));
			Vector temp=diff(*v3,(multiplyby_scalar(newV1,res1)));
			Vector newV3=diff(temp,(multiplyby_scalar(newV2,res2)));
			
			(*v2).x=newV2.x;
			(*v2).y=newV2.y;
			(*v2).z=newV2.z;
			
			(*v3).x=newV3.x;
			(*v3).y=newV3.y;
			(*v3).z=newV3.z;
		}
		else
		printf("v1 vektoru bos vektor.../n sifira bolunme hatasi...");
		
		a=1;	
		
	}
	else if (result==0)
	{
		a=0;
	}
	return a;
}
	 
char* vector2str(const Vector v) 
{
	char *str=malloc(sizeof(char)*15);
	sprintf(str,"(%.2f, %.2f, %.2f)",v.x,v.y,v.z);
	return str;
}
  
int main ()  
{ 
 Vector v1 = {1, 2, 2}, v2 = {-1, 0, 2}, v3 = {0, 0, 1}; 
 double k = 2; 
 
    printf("v1 = "); 
    print_vector(v1); 
    printf("v2 = "); 
    print_vector(v2); 
    printf("v3 = "); 
    print_vector(v3); 
   
    printf("v1 + v2 = "); 
 print_vector(sum(v1, v2)); 
  
  printf("v1 - v2 = "); 
 print_vector(diff(v1, v2)); 
  
 printf("k * v1  = "); 
 print_vector(multiplyby_scalar(v1, k)); 
  
 printf("v1 . v2 = %.2lf\n", dot_product(v1, v2)); 
  
 printf("v1 x v2 = "); 
print_vector(cross_product(v1, v2)); 
  
 printf("| v1 |  = %.2lf\n", norm(v1)); 
  
 if(is_unitvector(v1)) 
  printf("v1 is a unit vector.\n"); 
 else 
  printf("v1 is not unit vector.\n"); 
  
  printf("unit( v1 ) = "); 
 print_vector(unit(v1)); 
  
 printf("angle(v1, v2) = %.2lf\n", angle(v1, v2)); 
  
 printf("distance(v1, v2) = %.2lf\n", distance(v1, v2)); 
   
 if(are_linearly_independent(v1, v2, v3)) 
 printf("Vectors are linearly independent.\n"); 
 else 
 printf("Vectors are not linearly independent.\n"); 
  
 if(are_orthogonal(v1, v2, v3)) 
  printf("Vectors are orthogonal.\n"); 
 else 
  printf("Vectors are not orthogonal.\n"); 
   
 if(are_orthonormal(v1, v2, v3)) 
  printf("Vectors are orthonormal.\n"); 
 else 
  printf("Vectors are not orthonormal.\n");  
   
 printf("Projection of v1 onto v2 is = "); 
 print_vector(projection(v1, v2)); 
  
 printf("Orthogonal projection of v1 onto v2 is = "); 
 print_vector(orthogonal_projection(v1, v2)); 
 
  
  
 if(convert_2_orthogonal_basis(&v1, &v2, &v3)){ 
  printf("Orthogonalization of vectors:\n"); 
  printf("v1 = "); 
     print_vector(v1); 
     printf("v2 = "); 
     print_vector(v2); 
     printf("v3 = "); 
     print_vector(v3); 
 }      
    puts(vector2str(v1));  
    return 0;   
} 
