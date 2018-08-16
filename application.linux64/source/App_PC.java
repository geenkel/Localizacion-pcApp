import processing.core.*; 
import processing.data.*; 
import processing.event.*; 
import processing.opengl.*; 

import processing.serial.*; 
import java.math.*; 

import java.util.HashMap; 
import java.util.ArrayList; 
import java.io.File; 
import java.io.BufferedReader; 
import java.io.PrintWriter; 
import java.io.InputStream; 
import java.io.OutputStream; 
import java.io.IOException; 

public class App_PC extends PApplet {





Serial myPort;    // The serial port
PFont myFont;     // The display font
String inString;  // Input string from serial port
int lf = 10;      // ASCII linefeed 

//Distancias entre anclas
int d_ancla01=300;
int d_ancla02=384;
int d_ancla12=240;

int anchorPositioningFlag = 0;
int drawAnchorFlag = 0;

int [] center = {0,0};
//Coordenadas Ancla1 y Ancla2

int[] ancla1_ini = {0,0};
int[] ancla2_ini = {0,0};
int rotation = 0;
int[] ancla1 = {0,0};
int[] ancla2 = {0,0};

vec3d[] anchorArray = new vec3d[4];

Trilateration trilateration = new Trilateration();

int n_tags = 0;
ArrayList<Tag> tags = new ArrayList();
 
 //circulo de pruebas
int x = -200;
int y;
int r = 200;

int pos_ready = 0;
int errors = 0;
double error = 0;
vec3d T1 = new vec3d();
vec3d T2 = new vec3d();
int[] distances={0,0,0,0};

//*************************************************
public void setup() {
  
   
  center[0] = width/2;
  center[1] = height/2;
  
  myFont = loadFont("ArialMS-18.vlw"); 
  textFont(myFont, 18); 
  printArray(Serial.list());  
  
  //Encuentra las coordenadas de las anclas
  //anchor_positioning(d_ancla01,d_ancla02,d_ancla12);
  
  //Dibuja Mapa
  background(255);
  noStroke();
  fill(0);
  
  //anchorArray[0] = new vec3d(0,0,0);
  //anchorArray[1] = new vec3d(0,300,0);
  //anchorArray[2] = new vec3d(240,300,0);
  //anchorArray[3] = new vec3d(0,0,0);
  
  //draw_anchor();
  //draw_tag();
  
  delay(1000);
  //myPort = new Serial(this, Serial.list()[0], 115200); 
  myPort = new Serial(this, "/dev/ttyUSB0", 115200);
  myPort.bufferUntil(lf);  
      
  //vec3d tag= new vec3d(100,150,0);
  //tags.add(new Tag(4,(int)tag.x,(int)tag.y,0));
  print("Esperando posicionamiento de anclas");
} 

  

public void draw()
{
  if(pos_ready == 1)
  {
    pos_ready = 0;
    draw_tag();
  }
  if(drawAnchorFlag == 1)
  {
    drawAnchorFlag = 0;
    draw_anchor();
  }
}

public void serialEvent(Serial p) 
{
  inString = p.readString();
  print("string: "+inString);
  String[] distStr = inString.split(",");
  if(distStr.length>=3)
    distStr[3]=distStr[3].substring(0,distStr[3].length()-1);
  
  if(distStr.length>=4)
  {
    
    int[] idistances={0,0,0,0};
    vec3d solution= new vec3d();
    vec3d tag= new vec3d(x,y,0);
    
    //ID
    int id = Integer.parseInt(distStr[0]);
    println("ID="+id);
    
    //Distances
    idistances[0] = Integer.parseInt(distStr[1]);
    idistances[1] = Integer.parseInt(distStr[2]);
    idistances[2] = Integer.parseInt(distStr[3]);
    //println("d0 "+idistances[0]);
    //println("d1 "+idistances[1]);
    //println("d2 "+idistances[2]);
    
    if(id == 0)  //Posicionamiento de anclas
    {
      anchor_positioning(idistances[0],idistances[1],idistances[2]);
      drawAnchorFlag=1;
    }
    else if(id>2 && anchorPositioningFlag == 1)
    {
      boolean solution_found = false;
      for(int i = 0; i < 10 ; i++){
        trilateration.GetLocation2D(solution, 0 ,anchorArray, idistances);
        if(solution.x > 0.1f && solution.y > 0.1f){
          solution_found = true;
          break;
        }
        idistances[0] += 10;
        idistances[1] += 10;
        idistances[2] += 10;
      }  
      if(solution_found)
      {
        
          //Obtiene indice
        int index = -1;
        for(int i = 0; i<tags.size() ; i++)
        {
          if(tags.get(i).id == id)
          {
            index = i;
          }
        }
        if(index == -1){
          tags.add(new Tag(id, (int)solution.x, (int)solution.y, 200));
          index = tags.size() -1;
        }
        //tags.get(index).x = (int)solution.x;
        //tags.get(index).y = (int)solution.y;
        tags.get(index).locate((int)solution.x,(int)solution.y);
        tags.get(index).fcolor=200;
        println("("+tags.get(index).x+","+tags.get(index).y+")");
        ////arc((int)solution.x, (int)solution.y, 5, 5, 0, TWO_PI);
        pos_ready = 1;
      }
    }
    //idistances[0] = (200);
    //idistances[1] = (200);
    //idistances[2] = (360);
    
    //print("d1  "+idistances[0]+"\n");
    //print("d2  "+idistances[1]+"\n");
    //print("d3  "+idistances[2]+"\n");
    //print("d4  "+idistances[3]+"\n");
  }
  //print("\n Distancias ",d1,d2,d3);
}

public void keyPressed()
{
 if(key == 'R' || key == 'r')
 {
   rotate_map();
   draw_anchor();
 }
 if(key == 'M' || key == 'm')
 {
   mirror('x');
   draw_anchor();
 }
 if(key == 'N' || key == 'n')
 {
   mirror('y');
   draw_anchor();
 }
}

public void mirror(char axis)
{
  if(axis == 'x')
  {
    anchorArray[1].x = -anchorArray[1].x;
    anchorArray[2].x = -anchorArray[2].x;
    ancla1_ini[0] = -ancla1_ini[0];
    ancla2_ini[0] = -ancla2_ini[0];
  }
  else
  {
    anchorArray[1].y = -anchorArray[1].y;
    anchorArray[2].y = -anchorArray[2].y;
    ancla1_ini[1] = -ancla1_ini[1];
    ancla2_ini[1] = -ancla2_ini[1];
  }
}

public void rotate_map()
{
  rotation = rotation +10;
  if(rotation == 360)
    rotation=0;
  double theta = rotation * PI / 180;
  double x = ancla1_ini[0] * Math.cos(theta) - ancla1_ini[1] * Math.sin(theta);
  double y = ancla1_ini[0] * Math.sin(theta) + ancla1_ini[1] * Math.cos(theta);
  anchorArray[1].x = (int)(Math.round(x));
  anchorArray[1].y = (int)(Math.round(y));
  
  x = ancla2_ini[0] * Math.cos(theta) - ancla2_ini[1] * Math.sin(theta);
  y = ancla2_ini[0] * Math.sin(theta) + ancla2_ini[1] * Math.cos(theta);
  anchorArray[2].x = (int)(Math.round(x));
  anchorArray[2].y = (int)(Math.round(y));
}

public void draw_anchor()
{
  background(255);
  noStroke();
  fill(0);
  arc(center[0], center[1], 20, 20, 0, TWO_PI);
  text("A0 ", width/2+20, height/2+20);
  arc((int)anchorArray[1].x+center[0], (int)anchorArray[1].y+center[1], 20, 20, 0, TWO_PI);
  text("A1 ", (int)anchorArray[1].x+center[0]+20, (int)anchorArray[1].y+center[1]+20);
  arc((int)anchorArray[2].x+center[0], (int)anchorArray[2].y+center[1], 20, 20, 0, TWO_PI);
  text("A2 ", (int)anchorArray[2].x+center[0]+20, (int)anchorArray[2].y+center[1]+20);
  
  stroke(0);
  line(0, center[1] , width, center[1]);
  line(center[0],0,center[0],height);
}

public void draw_tag()
{
  draw_anchor();
  for(int i = 0; i < tags.size(); i++)
  {
    tags.get(i).draw(center[0],center[1]);
    text("T"+tags.get(i).id, tags.get(i).x+center[0]+20, tags.get(i).y+center[1]+20);
  }
}

class Anchor
{
  int id;
  int x;
  int y;
  public void draw(int fcolor, int cx, int cy)
  {  
    fill(100,100,100);
    arc(cx+x, cy+y, 20, 20, 0, TWO_PI);
  }
}


class Tag
{
  
  int id;
  int x;
  int y;
  int filter[][];
  int filter_count;
  boolean filter_init;
  int fcolor;
  int filter_size;
  
  Tag(int id, int x, int y, int fcolor)
  {
    this.id=id;
    this.x = x;
    this.y = y;
    this.fcolor = fcolor;
    filter_size = 3;
    filter = new int[filter_size][2];
    for(int i = 0; i < filter_size ; i++)
    {
      filter[i][0] = x;
      filter[i][1] = y;
    }
    filter_count = 0;
  }
  
  public void locate(int _x, int _y)
  {
    filter[filter_count][0] = _x;
    filter[filter_count][1] = _y;
    x = 0;
    y = 0;
    for(int i = 0; i <filter_size; i++)
    {
      x += filter[i][0];
      y += filter[i][1];
    }
    x /= filter_size;
    y /= filter_size;
    filter_count++;
    if(filter_count ==filter_size)
      filter_count = 0;
  }
  
  public void draw( int cx, int cy)
  {  
    fill(fcolor);
    arc(cx+x, cy+y, 10, 10, 0, TWO_PI);
  }
}

public void anchor_positioning(double d01, double d02, double d12)
{
  
  double r1 = d02;
  double r2 = d12;
  
  //Centro de la ancla 0
  double a1 = 0;
  double b1 = 0;
  
  //Centro de la ancla 1
  double a2 = d01;
  double b2 = 0;

  double D =-2*a1;
  double D2 =-2*a2;
  double E =-2*b1;
  double E2 =-2*b2;
  double F =a1*a1+b1*b1-r1*r1;
  double F2 =a2*a2+b2*b2-r2*r2;
  double G =(E2-E)/(D-D2);
  double H =(F2-F)/(D-D2);
  
  double A = G*G + 1;
  double B = G*D + 2*G*H + E;
  double C = H*H + H*D + F;
  double y1;
  double y2;
  
  y1 = ( -B + Math.sqrt(B*B-4*A*C) ) / ( 2*A );
  y2 = ( -B - Math.sqrt(B*B-4*A*C) ) / ( 2*A );
  
  A = 1;
  B = D2;
  C = y1*y1 + E2*y1 + F2;
  
  print("Las raíces son ",y1,y2);
  print("\n discriminante " , B*B-4*A*C);
  double x1 = ( -B + Math.sqrt(B*B-4*A*C) ) / ( 2*A );
  double x2 = ( -B - Math.sqrt(B*B-4*A*C) ) / ( 2*A );
  C = y2*y2 + E2*y2 + F2;
  print("\n discriminante " , B*B-4*A*C);
  double x3 = ( -B + Math.sqrt(B*B-4*A*C) ) / ( 2*A );
  double x4 = ( -B - Math.sqrt(B*B-4*A*C) ) / ( 2*A );
  
  ancla1_ini[0]= (int)a2;
  ancla1_ini[1]= (int)b2;
  ancla2_ini[0]= (int)x2;
  ancla2_ini[1]= (int)y1;
  
  ancla1[0]= (int)a2;
  ancla1[1]= (int)b2;
  ancla2[0]= (int)x2;
  ancla2[1]= (int)y1;
  
  
  //print("\n Las x son " , x1 , x2 , x3 , x4);
  
  //print("\n Las coordenadas son: ", x2 , y1+"\n");
  
  anchorArray[0] = new vec3d(0,0,0);
  anchorArray[1] = new vec3d(ancla1_ini[0],ancla1_ini[1],0);
  anchorArray[2] = new vec3d(ancla2_ini[0],ancla2_ini[1],0);
  anchorArray[3] = new vec3d(0,0,0);
  
  print("A1 = "+"("+ancla1[0]+","+ancla1[1]+")");
  print("A2 = "+"("+ancla2[0]+","+ancla2[1]+")");
  
  anchorPositioningFlag = 1;
}

//*******************************************************


//********************************************************

class vec3d
{
  double x;
  double y;
  double z;
  vec3d(){};
  vec3d(double x, double y, double z)
  {
    this.x=x;
    this.y=y;
    this.z=z;
  }
  vec3d(vec3d v)
  {
    this.x = v.x;
    this.y = v.y;
    this.z = v.z;
  }
  public void vprint()
  {
    //print("("+x+","+y+","+z+")\n");
  }
  public void copy(vec3d v)
  {
    x = v.x;
    y = v.y;
    z = v.z;
  }
}

class Trilateration
{
   public vec3d vdiff(vec3d vector1, vec3d vector2)
  {
    vec3d v= new vec3d();
    v.x = vector1.x - vector2.x;
    v.y = vector1.y - vector2.y;
    v.z = vector1.z - vector2.z;
    return v;
  }
  public vec3d vsum(vec3d vector1, vec3d vector2)
  {
    vec3d v= new vec3d();
    v.x = vector1.x + vector2.x;
    v.y = vector1.y + vector2.y;
    v.z = vector1.z + vector2.z;
    return v;
  }
  public vec3d vmul(vec3d vector, double n)
  {
    vec3d v= new vec3d();
    v.x = vector.x * n;
    v.y = vector.y * n;
    v.z = vector.z * n;
    return v;
  }
  public vec3d vdiv(vec3d vector, double n)
  {
    vec3d v = new vec3d();
    v.x = vector.x / n;
    v.y = vector.y / n;
    v.z = vector.z / n;
    return v;
  }
  public double vdist(vec3d v1, vec3d v2)
  {
    double xd = v1.x - v2.x;
    double yd = v1.y - v2.y;
    double zd = v1.z - v2.z;
    return Math.sqrt(xd * xd + yd * yd + zd * zd);
  }
  public double vnorm(vec3d vector)
  {
    return Math.sqrt(vector.x * vector.x + vector.y * vector.y + vector.z * vector.z);
  }
  public double dot(vec3d vector1, vec3d vector2)
  {
    return vector1.x * vector2.x + vector1.y * vector2.y + vector1.z * vector2.z;
  }
  public vec3d cross(vec3d vector1, vec3d vector2)
  {
    vec3d v = new vec3d();
    v.x = vector1.y * vector2.z - vector1.z * vector2.y;
    v.y = vector1.z * vector2.x - vector1.x * vector2.z;
    v.z = vector1.x * vector2.y - vector1.y * vector2.x;
    return v;
  }
  
  public double gdoprate(vec3d tag, vec3d p1, vec3d p2, vec3d p3)
  {
    vec3d ex;
    vec3d t1;
    vec3d t2;
    vec3d t3;
    
    double h, gdop1, gdop2, gdop3, result;
  
    ex = vdiff(p1, tag);
    h = vnorm(ex);
    t1 = vdiv(ex, h);
  
    ex = vdiff(p2, tag);
    h = vnorm(ex);
    t2 = vdiv(ex, h);
  
    ex = vdiff(p3, tag);
    h = vnorm(ex);
    t3 = vdiv(ex, h);
  
    gdop1 = Math.abs(dot(t1, t2));
    gdop2 = Math.abs(dot(t2, t3));
    gdop3 = Math.abs(dot(t3, t1));
  
    if (gdop1 < gdop2) result = gdop2; else result = gdop1;
    if (result < gdop3) result = gdop3;
  
    return result;
  }
  //**************************************************************************************************************************
  double mu1=0;
  double mu2=0;
  public int sphereline( vec3d p1,  vec3d p2,  vec3d sc, double r)
  {
      //println("P1="+p1.x+","+p1.y+" --- "+p2.x+","+p2.y+" --- "+sc.x+","+sc.y+" --- "+r);
     double a,b,c;
     double bb4ac;
     vec3d dp = new vec3d();
  
     dp.x = p2.x - p1.x;
     dp.y = p2.y - p1.y;
     dp.z = p2.z - p1.z;
  
     a = dp.x * dp.x + dp.y * dp.y + dp.z * dp.z;
  
     b = 2 * (dp.x * (p1.x - sc.x) + dp.y * (p1.y - sc.y) + dp.z * (p1.z - sc.z));
  
     c = sc.x * sc.x + sc.y * sc.y + sc.z * sc.z;
     c += p1.x * p1.x + p1.y * p1.y + p1.z * p1.z;
     c -= 2 * (sc.x * p1.x + sc.y * p1.y + sc.z * p1.z);
     c -= r * r;
  
     bb4ac = b * b - 4 * a * c;
  
     if (Math.abs(a) == 0 || bb4ac < 0) {
        mu1 = 0;
        mu2 = 0;
        return -1;
     }
  
     mu1 = (-b + Math.sqrt(bb4ac)) / (2 * a);
     mu2 = (-b - Math.sqrt(bb4ac)) / (2 * a);
  
     return 0;
  }
  //*************************************************************************************************************************************
  int    ERR_TRIL_CONCENTRIC           = -1;
  int    ERR_TRIL_COLINEAR_2SOLUTIONS   =   -2;
  int    ERR_TRIL_SQRTNEGNUMB         = -3;
  int    ERR_TRIL_NOINTERSECTION_SPHERE4  =    -4;
  int    ERR_TRIL_NEEDMORESPHERE        =  -5;
  int    TRIL_3SPHERES             = 3;
  int    TRIL_4SPHERES             = 4;
  public int trilateration(vec3d result1,
            vec3d result2,
            vec3d best_solution,
             vec3d p1,  double r1,
             vec3d p2,  double r2,
             vec3d p3,  double r3,
             vec3d p4,  double r4,
             double maxzero)
  {
    vec3d ex;
    vec3d ey;
    vec3d ez;
    vec3d t1;
    vec3d t2;
    vec3d t3;
    
    double  h, i, j, x, y, z, t;
    double mu;
    int result;
  
    /*********** FINDING TWO POINTS FROM THE FIRST THREE SPHERES **********/
  
    // if there are at least 2 concentric spheres within the first 3 spheres
    // then the calculation may not continue, drop it with error -1
  
    /* h = |p3 - p1|, ex = (p3 - p1) / |p3 - p1| */
    ex = vdiff(p3, p1); // vector p13
    h = vnorm(ex); // scalar p13
    if (h <= maxzero) {
      /* p1 and p3 are concentric, not good to obtain a precise intersection point */
      //printf("concentric13 return -1\n");
      return ERR_TRIL_CONCENTRIC;
    }
  
    /* h = |p3 - p2|, ex = (p3 - p2) / |p3 - p2| */
    ex = vdiff(p3, p2); // vector p23
    h = vnorm(ex); // scalar p23
    if (h <= maxzero) {
      /* p2 and p3 are concentric, not good to obtain a precise intersection point */
      //printf("concentric23 return -1\n");
      return ERR_TRIL_CONCENTRIC;
    }
  
    /* h = |p2 - p1|, ex = (p2 - p1) / |p2 - p1| */
    ex = vdiff(p2, p1); // vector p12
    h = vnorm(ex); // scalar p12
    if (h <= maxzero) {
      /* p1 and p2 are concentric, not good to obtain a precise intersection point */
      //printf("concentric12 return -1\n");
      return ERR_TRIL_CONCENTRIC;
    }
    ex = vdiv(ex, h); // unit vector ex with respect to p1 (new coordinate system)
  
    /* t1 = p3 - p1, t2 = ex (ex . (p3 - p1)) */
    t1 = vdiff(p3, p1); // vector p13
    i = dot(ex, t1); // the scalar of t1 on the ex direction
    t2 = vmul(ex, i); // colinear vector to p13 with the length of i
  
    /* ey = (t1 - t2), t = |t1 - t2| */
    ey = vdiff(t1, t2); // vector t21 perpendicular to t1
    t = vnorm(ey); // scalar t21
    if (t > maxzero) {
      /* ey = (t1 - t2) / |t1 - t2| */
      ey = vdiv(ey, t); // unit vector ey with respect to p1 (new coordinate system)
  
      /* j = ey . (p3 - p1) */
      j = dot(ey, t1); // scalar t1 on the ey direction
    } else
      j = 0.0f;
      
    //print("Dentro de Trilateracion= "+ex.x+","+ex.y+"--"+h+"--"+"t1="+t1.x+","+t1.y+"--"+"t2="+t2.x+","+t2.y+"---i="+i+"---t="+t+"---ey="+ey.x+","+ey.y+"\n");
    /* Note: t <= maxzero implies j = 0.0. */
    if (Math.abs(j) <= maxzero) {
      //println("J="+j);
      /* Is point p1 + (r1 along the axis) the intersection? */
      t2 = vsum(p1, vmul(ex, r1));
      if (Math.abs(vnorm(vdiff(p2, t2)) - r2) <= maxzero &&
          Math.abs(vnorm(vdiff(p3, t2)) - r3) <= maxzero) {
        /* Yes, t2 is the only intersection point. */
        result1.copy(t2);
        result2.copy(t2);
        return TRIL_3SPHERES;
      }
  
      /* Is point p1 - (r1 along the axis) the intersection? */
      t2 = vsum(p1, vmul(ex, -r1));
      if (Math.abs(vnorm(vdiff(p2, t2)) - r2) <= maxzero &&
          Math.abs(vnorm(vdiff(p3, t2)) - r3) <= maxzero) {
        /* Yes, t2 is the only intersection point. */
        result1.copy(t2);
        result2.copy(t2);
        return TRIL_3SPHERES;
      }
      /* p1, p2 and p3 are colinear with more than one solution */
      return ERR_TRIL_COLINEAR_2SOLUTIONS;
    }
  
    /* ez = ex x ey */
    ez = cross(ex, ey); // unit vector ez with respect to p1 (new coordinate system)
    //println("ez="+ez.x+","+ez.y);
  
    x = (r1*r1 - r2*r2) / (2*h) + h / 2;
    y = (r1*r1 - r3*r3 + i*i) / (2*j) + j / 2 - x * i / j;
    z = r1*r1 - x*x - y*y;
    if (z < -maxzero) {
      /* The solution is invalid, square root of negative number */
      return ERR_TRIL_SQRTNEGNUMB;
    } else
    if (z > 0.0f)
      z = Math.sqrt(z);
    else
      z = 0.0f;
  
    /* t2 = p1 + x ex + y ey */
    t2 = vsum(p1, vmul(ex, x));
    t2 = vsum(t2, vmul(ey, y));
  
    /* result1 = p1 + x ex + y ey + z ez */
    result1.copy(vsum(t2, vmul(ez, z)));
  
    /* result1 = p1 + x ex + y ey - z ez */
    result2.copy(vsum(t2, vmul(ez, -z)));
  
    /*********** END OF FINDING TWO POINTS FROM THE FIRST THREE SPHERES **********/
    /********* RESULT1 AND RESULT2 ARE SOLUTIONS, OTHERWISE RETURN ERROR *********/
    //println("Result1="+result1.x+"  "+result1.y);
    //println("Result2="+result2.x+"  "+result2.y);
  
    /************* FINDING ONE SOLUTION BY INTRODUCING ONE MORE SPHERE ***********/
  
    // check for concentricness of sphere 4 to sphere 1, 2 and 3
    // if it is concentric to one of them, then sphere 4 cannot be used
    // to determine the best solution and return -1
  
    /* h = |p4 - p1|, ex = (p4 - p1) / |p4 - p1| */
    ex = vdiff(p4, p1); // vector p14
    h = vnorm(ex); // scalar p14
    if (h <= maxzero) {
      /* p1 and p4 are concentric, not good to obtain a precise intersection point */
      //println("h <= maxzero");
      return TRIL_3SPHERES;
    }
    /* h = |p4 - p2|, ex = (p4 - p2) / |p4 - p2| */
    ex = vdiff(p4, p2); // vector p24
    h = vnorm(ex); // scalar p24
    if (h <= maxzero) {
      /* p2 and p4 are concentric, not good to obtain a precise intersection point */
      //printf("concentric24 return 0\n");
      return TRIL_3SPHERES;
    }
    /* h = |p4 - p3|, ex = (p4 - p3) / |p4 - p3| */
    ex = vdiff(p4, p3); // vector p34
    h = vnorm(ex); // scalar p34
    if (h <= maxzero) {
      /* p3 and p4 are concentric, not good to obtain a precise intersection point */
      //printf("concentric34 return 0\n");
      return TRIL_3SPHERES;
    }
  
    // if sphere 4 is not concentric to any sphere, then best solution can be obtained
    /* find i as the distance of result1 to p4 */
    t3 = vdiff(result1, p4);
    i = vnorm(t3);
    /* find h as the distance of result2 to p4 */
    t3 = vdiff(result2, p4);
    h = vnorm(t3);
  
    /* pick the result1 as the nearest point to the center of sphere 4 */
    if (i > h) {
      best_solution.copy(result1);
      result1.copy(result2);
      result2.copy(best_solution);
    }
  
    int count4 = 0;
    double rr4 = r4;
    result = 1;
    /* intersect result1-result2 vector with sphere 4 */
    while(result!=0 && count4 < 10)
    {
      result=sphereline(result1, result2, p4, rr4);
      rr4+=0.1f;
      count4++;
    }
  
    if (result != 0) {
  
      /* No intersection between sphere 4 and the line with the gradient of result1-result2! */
      best_solution.copy(result1); // result1 is the closer solution to sphere 4
      //return ERR_TRIL_NOINTERSECTION_SPHERE4;
  
    } else {
  
      if (mu1 < 0 && mu2 < 0) {
  
        /* if both mu1 and mu2 are less than 0 */
        /* result1-result2 line segment is outside sphere 4 with no intersection */
        if (Math.abs(mu1) <= Math.abs(mu2)) mu = mu1; else mu = mu2;
        /* h = |result2 - result1|, ex = (result2 - result1) / |result2 - result1| */
        ex = vdiff(result2, result1); // vector result1-result2
        h = vnorm(ex); // scalar result1-result2
        ex = vdiv(ex, h); // unit vector ex with respect to result1 (new coordinate system)
        /* 50-50 error correction for mu */
        mu = 0.5f*mu;
        /* t2 points to the intersection */
        t2 = vmul(ex, mu*h);
        t2 = vsum(result1, t2);
        /* the best solution = t2 */
        best_solution.copy(t2);
  
      } else if ((mu1 < 0 && mu2 > 1) || (mu2 < 0 && mu1 > 1)) {
  
        /* if mu1 is less than zero and mu2 is greater than 1, or the other way around */
        /* result1-result2 line segment is inside sphere 4 with no intersection */
        if (mu1 > mu2) mu = mu1; else mu = mu2;
        /* h = |result2 - result1|, ex = (result2 - result1) / |result2 - result1| */
        ex = vdiff(result2, result1); // vector result1-result2
        h = vnorm(ex); // scalar result1-result2
        ex = vdiv(ex, h); // unit vector ex with respect to result1 (new coordinate system)
        /* t2 points to the intersection */
        t2 = vmul(ex, mu*h);
        t2 = vsum(result1, t2);
        /* vector t2-result2 with 50-50 error correction on the length of t3 */
        t3 = vmul(vdiff(result2, t2),0.5f);
        /* the best solution = t2 + t3 */
        best_solution.copy(vsum(t2, t3));
  
      } else if (((mu1 > 0 && mu1 < 1) && (mu2 < 0 || mu2 > 1))
          || ((mu2 > 0 && mu2 < 1) && (mu1 < 0 || mu1 > 1))) {
  
        /* if one mu is between 0 to 1 and the other is not */
        /* result1-result2 line segment intersects sphere 4 at one point */
        if (mu1 >= 0 && mu1 <= 1) mu = mu1; else mu = mu2;
        /* add or subtract with 0.5*mu to distribute error equally onto every sphere */
        if (mu <= 0.5f) mu-=0.5f*mu; else mu-=0.5f*(1-mu);
        /* h = |result2 - result1|, ex = (result2 - result1) / |result2 - result1| */
        ex = vdiff(result2, result1); // vector result1-result2
        h = vnorm(ex); // scalar result1-result2
        ex = vdiv(ex, h); // unit vector ex with respect to result1 (new coordinate system)
        /* t2 points to the intersection */
        t2 = vmul(ex, mu*h);
        t2 = vsum(result1, t2);
        /* the best solution = t2 */
        best_solution.copy(t2);
  
      } else if (mu1 == mu2) {
  
        /* if both mu1 and mu2 are between 0 and 1, and mu1 = mu2 */
        /* result1-result2 line segment is tangential to sphere 4 at one point */
        mu = mu1;
        /* add or subtract with 0.5*mu to distribute error equally onto every sphere */
        if (mu <= 0.25f) mu-=0.5f*mu;
        else if (mu <=0.5f) mu-=0.5f*(0.5f-mu);
        else if (mu <=0.75f) mu-=0.5f*(mu-0.5f);
        else mu-=0.5f*(1-mu);
        /* h = |result2 - result1|, ex = (result2 - result1) / |result2 - result1| */
        ex = vdiff(result2, result1); // vector result1-result2
        h = vnorm(ex); // scalar result1-result2
        ex = vdiv(ex, h); // unit vector ex with respect to result1 (new coordinate system)
        /* t2 points to the intersection */
        t2 = vmul(ex, mu*h);
        t2 = vsum(result1, t2);
        /* the best solution = t2 */
        best_solution.copy(t2);
  
      } else {
  
        /* if both mu1 and mu2 are between 0 and 1 */
        /* result1-result2 line segment intersects sphere 4 at two points */
  
        //return ERR_TRIL_NEEDMORESPHERE;
  
        mu = mu1 + mu2;
        /* h = |result2 - result1|, ex = (result2 - result1) / |result2 - result1| */
        ex = vdiff(result2, result1); // vector result1-result2
        h = vnorm(ex); // scalar result1-result2
        ex = vdiv(ex, h); // unit vector ex with respect to result1 (new coordinate system)
        /* 50-50 error correction for mu */
        mu = 0.5f*mu;
        /* t2 points to the intersection */
        t2 = vmul(ex, mu*h);
        t2 = vsum(result1, t2);
        /* the best solution = t2 */
        best_solution.copy(t2);
  
      }
  
    }
  
    return TRIL_4SPHERES;
  
    /******** END OF FINDING ONE SOLUTION BY INTRODUCING ONE MORE SPHERE *********/
  }
  //***********************************************************************************************************************************************************
  int nosolution_count=0;
  double best_3derror=0;
  double best_gdoprate=0;
  int combination=0;
  double MAXZERO = 0.001f;
  public int deca_3dlocate (   vec3d solution1,
                        vec3d solution2,
                        vec3d best_solution,
                        vec3d p1_, double r1,
                        vec3d p2_, double r2,
                        vec3d p3_, double r3,
                        vec3d p4_, double r4)
  {
    vec3d p1 = new vec3d(p1_);
    vec3d p2 = new vec3d(p2_);
    vec3d p3 = new vec3d(p3_);
    vec3d p4 = new vec3d(p4_);
    
    vec3d o1 = new vec3d();
    vec3d o2 = new vec3d();
    vec3d solution = new vec3d();
    vec3d ptemp = new vec3d();
    
    vec3d  solution_compare1 = new vec3d();
    vec3d solution_compare2 = new vec3d();
    
    double  /*error_3dcompare1, error_3dcompare2,*/ rtemp;
    double  gdoprate_compare1, gdoprate_compare2;
    double  ovr_r1, ovr_r2, ovr_r3, ovr_r4;
    int    overlook_count, combination_counter;
    int    trilateration_errcounter, trilateration_mode34;
    int    success, concentric, result;
  
    trilateration_errcounter = 0;
    trilateration_mode34 = 0;
  
    combination_counter = 4; /* four spheres combination */
  
    best_gdoprate = 1; /* put the worst gdoprate init */
    gdoprate_compare1 = 1; gdoprate_compare2 = 1;
    solution_compare1.x = 0; solution_compare1.y = 0; solution_compare1.z = 0;
      //error_3dcompare1 = 0;
  
    do {
      success = 0;
      concentric = 0;
      overlook_count = 0;
      ovr_r1 = r1; ovr_r2 = r2; ovr_r3 = r3; ovr_r4 = r4;
      //print("\nOVR = "+ovr_r1+"\t"+ovr_r2+"\t"+ovr_r3+"\t"+ovr_r4+"\n");
      //print("Overlook="+overlook_count+"\n");
      //print("Concentric="+concentric+"\n");
      do {
  
        result = trilateration(o1, o2, solution, p1, ovr_r1, p2, ovr_r2, p3, ovr_r3, p4, ovr_r4, MAXZERO);
        //print("Result="+result+"\n");
        //print("o1=");
        o1.vprint();
        //print("o2=");
        o2.vprint();
        //print("Solution=");
        solution.vprint();
        switch (result)
        {
          case 3: // 3 spheres are used to get the result
            trilateration_mode34 = TRIL_3SPHERES;
            success = 1;
            break;
  
          case 4: // 4 spheres are used to get the result
            trilateration_mode34 = TRIL_4SPHERES;
            success = 1;
            break;
  
          case -1:
            concentric = 1;
            break;
  
          default: // any other return value goes here
            ovr_r1 += 0.10f;
            ovr_r2 += 0.10f;
            ovr_r3 += 0.10f;
            ovr_r4 += 0.10f;
            overlook_count++;
            break;
        }
  
              //qDebug() << "while(!success)" << overlook_count << concentric << "result" << result;
  
      } while (success==0 && (overlook_count <= 5) && concentric==0);
  
  
      if (success!=0)
      {
        switch (result)
        {
        case 3:
          solution1.copy(o1);
          solution2.copy(o2);
          //println("Solutions=");
          //o1.vprint();
          //o2.vprint();
          nosolution_count = overlook_count;
          combination_counter = 0;
          break;
  
        case 4:      //TRIL_4SPHERES
          /* calculate the new gdop */
          gdoprate_compare1  = gdoprate(solution, p1, p2, p3);
  
          /* compare and swap with the better result */
          if (gdoprate_compare1 <= gdoprate_compare2) {
  
            solution1.copy(o1);
            solution2.copy(o2);
            best_solution.copy(solution);
            nosolution_count = overlook_count;
            best_3derror  = Math.sqrt((vnorm(vdiff(solution, p1))-r1)*(vnorm(vdiff(solution, p1))-r1) +
                      (vnorm(vdiff(solution, p2))-r2)*(vnorm(vdiff(solution, p2))-r2) +
                      (vnorm(vdiff(solution, p3))-r3)*(vnorm(vdiff(solution, p3))-r3) +
                      (vnorm(vdiff(solution, p4))-r4)*(vnorm(vdiff(solution, p4))-r4));
            best_gdoprate  = gdoprate_compare1;
  
            /* save the previous result */
            solution_compare2.copy(solution_compare1); 
                      //error_3dcompare2 = error_3dcompare1;
            gdoprate_compare2 = gdoprate_compare1;
  
            combination = 5 - combination_counter;
  
            ptemp.copy(p1); p1.copy(p2); p2.copy(p3); p3.copy(p4); p4.copy(ptemp);
            rtemp = r1; r1 = r2; r2 = r3; r3 = r4; r4 = rtemp;
            combination_counter--;
  
          }
          break;
  
        default:
          break;
        }
      }
      else
      {
              //trilateration_errcounter++;
              trilateration_errcounter = 4;
              combination_counter = 0;
      }
  
          //ptemp = p1; p1 = p2; p2 = p3; p3 = p4; p4 = ptemp;
          //rtemp = r1; r1 = r2; r2 = r3; r3 = r4; r4 = rtemp;
          //combination_counter--;
          //qDebug() << "while(combination_counter)" << combination_counter;
  
      } while (combination_counter!=0);
  
    // if it gives error for all 4 sphere combinations then no valid result is given
    // otherwise return the trilateration mode used
    if (trilateration_errcounter >= 4) return -1; else return trilateration_mode34;
  
  }
  //***********************************************************************************************************************************************************************
  
  public int GetLocation2D(vec3d best_solution, int use4thAnchor, vec3d[] anchorArray, int[] distanceArray)
  {    
    //print("Distancias="+distanceArray[0],distanceArray[1],distanceArray[2]);
  
    vec3d o1 = new vec3d();
    vec3d o2 = new vec3d();
    vec3d p1 = new vec3d();
    vec3d p2 = new vec3d();
    vec3d p3 = new vec3d();
    vec3d p4 = new vec3d();
    double  r1 = 0, r2 = 0, r3 = 0, r4 = 0;
    int    result;
    int     error ;
  
    vec3d  t3 = new vec3d();
    double  dist1, dist2;
  
    /* Anchors coordinate */
      p1.x = anchorArray[0].x;    p1.y = anchorArray[0].y;  p1.z = anchorArray[0].z;
      p2.x = anchorArray[1].x;    p2.y = anchorArray[1].y;  p2.z = anchorArray[1].z;
      p3.x = anchorArray[2].x;    p3.y = anchorArray[2].y;  p3.z = anchorArray[2].z;
      p4.x = anchorArray[0].x;    p4.y = anchorArray[0].y;  p4.z = anchorArray[0].z; //4th same as 1st - only 3 used for trilateration
    
      //print("\nAnchor 0="); p1.vprint();
      //print("Anchor 1="); p2.vprint();
      //print("Anchor 2="); p3.vprint();
      //print("Anchor 3="); p4.vprint();
    
      r1 = (double) distanceArray[0] ;
      r2 = (double) distanceArray[1] ;
      r3 = (double) distanceArray[2] ;
      r4 = (double) distanceArray[0] ;
  
      //print("Distancias = ",r1,",",r2,",",r3,",",r4);
  
      r4 = r1;
  
    /* get the best location using 3 or 4 spheres and keep it as know_best_location */
      result = deca_3dlocate (o1, o2, best_solution,p1, r1, p2, r2, p3, r3, p4, r1);
      error = nosolution_count;
      //cout << "\n GetLocation   " << result <<" , " << "  sol1: " << o1.x <<" , " << o1.y <<" , " << o1.z <<" , " << " sol2: " << o2.x <<" , " << o2.y <<" , " << o2.z<<endl;
      //print( "\n GetLocation   " , result ," , " , "  sol1: " , o1.x ," , " , o1.y ," , " , o1.z ," , " , " sol2: " , o2.x ," , " , o2.y ," , " , o2.z,"\n");
      //print("\nError ="+error+"\n");
      best_solution.vprint();
  
      if(result >= 0)
      {
          if (use4thAnchor == 1) //if have 4 ranging results, then use 4th anchor to pick solution closest to it
          {
                  double diff1, diff2;
                  /* find dist1 as the distance of o1 to known_best_location */
                  t3 = vdiff(o1, anchorArray[3]);
                  dist1 = vnorm(t3);
  
                  t3 = vdiff(o2, anchorArray[3]);
                  dist2 = vnorm(t3);
  
                  /* find the distance closest to received range measurement from 4th anchor */
                  diff1 = Math.abs(r4 - dist1);
                  diff2 = Math.abs(r4 - dist2);
  
                  /* pick the closest match to the 4th anchor range */
                  if (diff1 < diff2) best_solution.copy(o1); else best_solution.copy(o2);
          }
          else
          {
              //assume tag is below the anchors (1, 2, and 3)
              if(o1.z < p1.z) best_solution.copy(o1); else best_solution.copy(o2);
          }
      }
  
    if (result >= 0)
    {
      return result;
    }
  
    //return error
    return -1;
  } 
}
  public void settings() {  size(1200,1000); }
  static public void main(String[] passedArgs) {
    String[] appletArgs = new String[] { "App_PC" };
    if (passedArgs != null) {
      PApplet.main(concat(appletArgs, passedArgs));
    } else {
      PApplet.main(appletArgs);
    }
  }
}
