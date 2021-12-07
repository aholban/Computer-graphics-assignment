#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <glm/glm.hpp>
#include "SDL.h"
#include <DrawingWindow.h>
#include <Colour.h>
#include <CanvasTriangle.h>
#include <Utils.h>
#include <ModelTriangle.h>
#include <cmath>

using namespace std;
using namespace glm;

#define WIDTH 640
#define HEIGHT 480

DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);


struct face{
  int pointLocation[3];
  Colour mat;
  int textureLocation[3];
};


void update();
void handleEvent(SDL_Event event, vec3 &cameraLocation, mat3 &cameraOrientation, vec2 &moveCamera, bool &filled, bool &texture, int &ppm);
CanvasPoint createVertice(CanvasPoint from, CanvasPoint to, CanvasPoint yCoord);
void drawLine(CanvasPoint from, CanvasPoint to, Colour colour, float buffer[WIDTH][HEIGHT]);
void drawTriangle(CanvasTriangle triangle, Colour colour, float buffer[WIDTH][HEIGHT]);
void drawFilledTriangle(CanvasTriangle triangle, Colour colour, float buffer[WIDTH][HEIGHT]);
void fillTopTriangle(CanvasTriangle triangle, Colour colour, float buffer[WIDTH][HEIGHT]);
void fillBottomTriangle(CanvasTriangle triangle, Colour colour, float buffer[WIDTH][HEIGHT]);
void readPPMFile(string filename, vector<Colour> &textureMap, int &textureHeight, int &textureWidth);
void readMtlFile(string filename, vector<Colour> &palette, vector<Colour> &textureMap, int &textureHeight, int &textureWidth);
void readObjFile(vector<ModelTriangle> &triangles,   vector<Colour> &palette, vector<Colour> &textureMap, int &textureHeight, int &textureWidth);
vector<CanvasTriangle>  getCanvas(vector<ModelTriangle> triangles, vec3 cameraLocation, mat3 cameraOrientation, vec2 moveCamera, int textureHeight, int textureWidth);
void rasterizer(vector<CanvasTriangle> allTri, float buffer[WIDTH][HEIGHT]);
void savePPMfile(unsigned int height, unsigned int width, unsigned int colour_max, int &ppm);
void textureMapping(vector<CanvasTriangle>allTris, vector<Colour>textureMap, int textureHeight, int textureWidth, float buffer[WIDTH][HEIGHT]);

void update()
{

  // Function for performing animation (shifting artifacts or moving the camera)
}

void handleEvent(SDL_Event event, vec3 &cameraLocation, mat3 &cameraOrientation, vec2 &moveCamera,bool &filled, bool &texture, int &ppm){
  if(event.type == SDL_KEYDOWN) {
    if(event.key.keysym.sym == SDLK_LEFT) {

      moveCamera.x -= 10;
    }
    else if(event.key.keysym.sym == SDLK_RIGHT) {

      moveCamera.x += 10;
    }
    else if(event.key.keysym.sym == SDLK_UP) {

      moveCamera.y += 10;
    }
    else if(event.key.keysym.sym == SDLK_DOWN) {

      moveCamera.y -= 10;
    }
    else if(event.key.keysym.sym == SDLK_w)
    {
      //if(cameraLocation.z + 2 >= -46){
        // cout<<"Too close"<<endl;
      //}
      //else{
        cameraLocation.z += 2;
      //}

    }
    else if(event.key.keysym.sym == SDLK_s)
    {
      cameraLocation.z -= 2;
    }
    else if(event.key.keysym.sym == SDLK_p){
      savePPMfile(HEIGHT,WIDTH,255,ppm);
    }
    else if(event.key.keysym.sym == SDLK_t){
      if(!texture) texture = true;
      else texture = false;
    }
    else if(event.key.keysym.sym == SDLK_f){
      if(!filled) filled = true;
      else filled = false;
    }
    else if(event.key.keysym.sym == SDLK_z){
      float angle = 0.01;
      mat3 rotationMatrix(vec3(cos(angle), -sin(angle), 0.0), vec3(sin(angle), cos(angle), 0.0), vec3(0.0, 0.0, 1.0));
      cameraOrientation = cameraOrientation * rotationMatrix;
      //while(cameraLocation.z >= -44) cameraLocation.z -= 2;
    }
    else if(event.key.keysym.sym == SDLK_x){
      float angle = 0.01;
      mat3 rotationMatrix(vec3(1.0, 0.0, 0.0), vec3(0.0, cos(angle), -sin(angle)), vec3(0.0, sin(angle), cos(angle)));
      cameraOrientation = cameraOrientation * rotationMatrix;
      //while(cameraLocation.z >= -44) cameraLocation.z -= 2;
    }
    else if(event.key.keysym.sym == SDLK_y){
      float angle = 0.01;
      mat3 rotationMatrix(vec3(cos(angle), 0.0, sin(angle)), vec3(0.0, 1.0, 0.0), vec3(-sin(angle), 0.0, cos(angle)));
      cameraOrientation = cameraOrientation * rotationMatrix;
    //  while(cameraLocation.z >= -44) cameraLocation.z -= 2;
    }
  }
}

CanvasPoint createVertice(CanvasPoint v1, CanvasPoint v2, CanvasPoint v3){
  float t = (v2.y - v1.y)/(v3.y - v1.y);
  CanvasPoint vertice = CanvasPoint(v1.x + t*(v3.x - v1.x), v2.y);
  float t1 = (vertice.x - v1.x)/(v3.x - v1.x);
  float z = v1.depth + t1*(v3.depth - v1.depth);
  vertice.depth = z;    //calculate the depth of the vertice
  float t2 = (v2.texturePoint.y - v1.texturePoint.y)/(v3.texturePoint.y - v1.texturePoint.y);
  vertice.texturePoint = TexturePoint(v1.texturePoint.x + t2*(v3.texturePoint.x - v1.texturePoint.x), v2.texturePoint.y);     //calculate the respective texture point
  return vertice;
}


void drawLine(CanvasPoint from, CanvasPoint to, Colour colour, float buffer[WIDTH][HEIGHT])
{
  //window.clearPixels();
  float xDiff = to.x - from.x;
  float yDiff = to.y - from.y;
  float zDiff = to.depth - from.depth;
  float numberOfSteps = std::max(abs(xDiff), abs(yDiff));
  float xStepSize = xDiff/numberOfSteps;
  float yStepSize = yDiff/numberOfSteps;
  float zStepSize = zDiff/numberOfSteps;
  uint32_t colourPixel = (255<<24) + (int(colour.red)<<16) + (int( colour.green)<<8) + int(colour.blue);
  for(float i = 0.0; i<numberOfSteps; i++){
    float x = from.x + (xStepSize * i);
    float y = from.y + (yStepSize * i);

    float z = from.depth +(zStepSize * i);

    if(x>=0 && y>=0 && x<=WIDTH && y<=HEIGHT){
      if(buffer[(int)round(x)][(int)round(y)] >= 1/z ) {
        buffer[(int)round(x)][(int)round(y)] = 1/z;
        window.setPixelColour(round(x), round(y), colourPixel);
      }
    }
  }
}

void drawTriangle(CanvasTriangle triangle, Colour colour, float buffer[WIDTH][HEIGHT])
{
  drawLine(triangle.vertices[0], triangle.vertices[1], colour, buffer);
  drawLine(triangle.vertices[0], triangle.vertices[2], colour, buffer);
  drawLine(triangle.vertices[1], triangle.vertices[2], colour, buffer);
}

CanvasTriangle sortVertices(CanvasTriangle triangle)
{
  if(triangle.vertices[0].y > triangle.vertices[1].y) swap(triangle.vertices[0], triangle.vertices[1]);
  if(triangle.vertices[0].y > triangle.vertices[2].y) swap(triangle.vertices[0], triangle.vertices[2]);
  if(triangle.vertices[1].y > triangle.vertices[2].y) swap(triangle.vertices[1], triangle.vertices[2]);
  return triangle;
}

void fillTopTriangle(CanvasTriangle triangle, Colour colour, float buffer[WIDTH][HEIGHT])
{
  float invslope1 = (triangle.vertices[1].x - triangle.vertices[0].x) / (triangle.vertices[1].y - triangle.vertices[0].y);
  float invslope2 = (triangle.vertices[2].x - triangle.vertices[0].x) / (triangle.vertices[2].y - triangle.vertices[0].y);

  float invslopeZ1 = (triangle.vertices[1].depth - triangle.vertices[0].depth) / (triangle.vertices[1].y - triangle.vertices[0].y);
  float invslopeZ2 = (triangle.vertices[2].depth - triangle.vertices[0].depth) / (triangle.vertices[2].y - triangle.vertices[0].y);

  float x1 = triangle.vertices[0].x;
  float x2 = triangle.vertices[0].x;

  float z1 = triangle.vertices[0].depth;
  float z2 = triangle.vertices[0].depth;

  for(int y = triangle.vertices[0].y; y<=triangle.vertices[1].y; y++)
  {
    CanvasPoint v1 = CanvasPoint(x1, y, z1);
    CanvasPoint v2 = CanvasPoint(x2, y, z2);
    drawLine(v1, v2, colour, buffer);
    x1 += invslope1;
    x2 += invslope2;
    z1 += invslopeZ1;
    z2 += invslopeZ2;
  }
}

void fillBottomTriangle(CanvasTriangle triangle, Colour colour, float buffer[WIDTH][HEIGHT])
{
  float invslope1 = (triangle.vertices[2].x - triangle.vertices[0].x) / (triangle.vertices[2].y - triangle.vertices[0].y);
  float invslope2 = (triangle.vertices[2].x - triangle.vertices[1].x) / (triangle.vertices[2].y - triangle.vertices[1].y);

  float invslopeZ1 = (triangle.vertices[2].depth - triangle.vertices[0].depth) / (triangle.vertices[2].y - triangle.vertices[0].y);
  float invslopeZ2 = (triangle.vertices[2].depth - triangle.vertices[1].depth) / (triangle.vertices[2].y - triangle.vertices[1].y);

  float x1 = triangle.vertices[2].x;
  float x2 = triangle.vertices[2].x;

  float z1 = triangle.vertices[2].depth;
  float z2 = triangle.vertices[2].depth;

  for(int y = triangle.vertices[2].y; y>triangle.vertices[0].y; y--)
  {
    CanvasPoint v1 = CanvasPoint(x1, y, z1);
    CanvasPoint v2 = CanvasPoint(x2, y, z2);
    drawLine(v1, v2, colour, buffer);
    x1 -= invslope1;
    x2 -= invslope2;
    z1 -= invslopeZ1;
    z2 -= invslopeZ2;
  }
}

void drawFilledTriangle(CanvasTriangle triangle, Colour colour, float buffer[WIDTH][HEIGHT])
{
  triangle = sortVertices(triangle);
  CanvasPoint vertice = createVertice(triangle.vertices[0],triangle.vertices[1],triangle.vertices[2]);
  CanvasTriangle triangle1 = CanvasTriangle(triangle.vertices[0],vertice,triangle.vertices[1]);
  CanvasTriangle triangle2 = CanvasTriangle(vertice,triangle.vertices[1],triangle.vertices[2]);
  fillTopTriangle(triangle1, colour, buffer);
  fillBottomTriangle(triangle2, colour, buffer);
}

void readPPMFile(string filename, vector<Colour> &textureMap, int &textureHeight, int &textureWidth){
  string line;

  // filename can come with a trailing \r, to remove this
  // or other newline chars we filter them out.
  string filteredFilename;
  for ( unsigned long int i = 0 ; i < filename.length() ; i++ )
  {
      if(!(filename[i] == '\n' || filename[i] == '\r')){
        filteredFilename.push_back(filename[i]);
      }
  }
  ifstream texture (filteredFilename,  std::ifstream::binary);

  // skip first 2 lines, they are comments
  getline(texture, line,'\n');
  //cout << line << endl;
  getline(texture, line,'\n');
  //cout << line << endl;

  // 3rd line will be the Height and Width of the PPM
  getline(texture, line,'\n');
  //cout << line << endl;
  string *words = split(line,' ');
  textureHeight = stoi(words[0]);
  textureWidth = stoi(words[1]);

  getline(texture, line,'\n');
  //cout << line << endl;
  //int colorRange = stoi(line);

  // Set of ints for the colours
  unsigned int red=0;
  unsigned int green=0;
  unsigned int blue=0;

  // Uint8 to hold the Char (8bits) directly from the file
  char charC[3];
  //track what pixel we are up to
  int x = 0;
  int y = 0;
  while(!texture.eof()){
    texture.read(charC, 3);
    // convert Char to Int
    red = int(charC[0]);
    green = int(charC[1]);
    blue = int(charC[2]);

    textureMap.push_back(Colour(red, green, blue));

    //update window
    //cout<<"At "<<x<<" and "<<y<<" colours: "<<red<<" "<<green<<" "<<blue<<endl;
  //  uint32_t colourPixel = (255<<24) + (red<<16) + (green<<8) + blue;
  //  window.setPixelColour(x, y, colourPixel);


    // iterate the pixel
    if(x == textureWidth - 1){
      // ond of line, should have a whitespace char to remove
      x = 0;
      y++;
    }else{
      x++;
    }
  }
  // printf("The end\n");
  texture.close();

}

void readMtlFile(const string filename, vector<Colour> &palette, vector<Colour> &textureMap, int &textureHeight, int &textureWidth)
{

string filteredFilename;
for ( unsigned long int i = 0 ; i < filename.length() ; i++ )
{
    if(!(filename[i] == '\n' || filename[i] == '\r')){
      filteredFilename.push_back(filename[i]);
    }
}
  ifstream mtlfile(filteredFilename, std::ifstream::in);
  string line, kd;
//  int r, g, b;
  Colour holder;

  while(!mtlfile.eof()){
    getline(mtlfile,line);
    string *words = split(line,' ');
    if(words[0]=="newmtl") {
      holder.name = words[1];
    }
    if(words[0]=="Kd"){

      holder.red = stof(words[1]) * 255.0f;
      holder.green = stof(words[2]) * 255.0f;
      holder.blue = stof(words[3]) * 255.0f;

      //cout << holder.name << " "<< holder.red << " "<< holder.green << " "<< holder.blue << " " <<endl;

      palette.push_back(holder);
    }
    if(words[0]=="map_Kd"){
      readPPMFile(words[1], textureMap, textureHeight, textureWidth);

    }
  }
  mtlfile.close();
}


void readObjFile(vector<ModelTriangle> &triangles,   vector<Colour> &palette, vector<Colour> &textureMap, int &textureHeight, int &textureWidth)
{
  ifstream file("logo.obj");
  string line;
  Colour material = Colour(255,255,255);
  vector<vec3> points;
  vector<vec2> texturePoints;
  vector<face> faces;
  int loc = 0;


  while(!file.eof()){
    getline(file,line);

    string *words = split(line,' ');  //split line into 'word' tokens
    //see what the line is
    if (words[0] == "o"){
      material = Colour();

      //cout << "New object being created" << endl;

    }else if(words[0]=="usemtl") {  //if a material declaration
      bool found = false;
      loc = 0;
      while(!found){  //loop through palette to see if the material is known to use
        if (loc >= (signed int)palette.size()){
          found = true;
          loc = -1; //incase it's not, don't just break

          //cout << "No mat found of name: " << words[1] << endl;
        }else{

          Colour toCheck = palette.at(loc);

          if (toCheck.name == words[1] ){ //if we find it, store it as current material
            found = true;
            material = toCheck;

            //cout << "new material: " << material.name << endl;
          } else {
            loc++;
          }

        }
      }

    } else if(words[0]=="v"){ //if a vector

      vec3 toBePlaced = vec3(stof(words[1]), stof(words[2]), stof(words[3]));
      points.push_back(toBePlaced);
      //cout<< words[1] << " " << words[2] << " " << words[3] << endl;

    }else if(words[0]=="vt"){  //if a texture vertice declaration store the coordinates in a vector

      vec2 toBePlaced = vec2(stof(words[1]), stof(words[2]));
      texturePoints.push_back(toBePlaced);


    }else if(words[0]=="f"){  //if a face declaration
      face newFace;
      string *number1 = split(words[1],'/');
      string *number2 = split(words[2],'/');
      string *number3 = split(words[3],'/');
      newFace.pointLocation[0] = (stoi(number1[0]));
      newFace.pointLocation[1] = (stoi(number2[0]));
      newFace.pointLocation[2] = (stoi(number3[0]));
      newFace.mat = material;

      //store the location in the texturePoints vector of the texture points for a face
      newFace.textureLocation[0] = (stoi(number1[1]));
      newFace.textureLocation[1] = (stoi(number2[1]));
      newFace.textureLocation[2] = (stoi(number3[1]));

      faces.push_back(newFace);

    //  cout << "Face at: " << newFace.pointLocation[0] << " "<< newFace.pointLocation[1] << " " << newFace.pointLocation[2] << ", with material: " << newFace.mat.name << endl;
  }else if(words[0]=="mtllib"){
      // cout<<"Got here"<<endl;
      readMtlFile(words[1], palette, textureMap, textureHeight, textureWidth);
    }
  }
  file.close();

  for (int i = 0; i < (signed int)faces.size(); i++){

    int xLoc = faces.at(i).pointLocation[0] ;
    int yLoc = faces.at(i).pointLocation[1] ;
    int zLoc = faces.at(i).pointLocation[2] ;
    Colour material = faces.at(i).mat;
    float scale = 0.05;
    vec3 scalar = vec3(scale, scale, 1);
    vec3 pointOne =  points.at(xLoc - 1) * scalar;
    vec3 pointTwo =  points.at(yLoc - 1) * scalar;
    vec3 pointThree =points.at(zLoc - 1) * scalar;

    //store the coordinates of the texture points for a face
    int texLoc0 = faces.at(i).textureLocation[0];
    int texLoc1 = faces.at(i).textureLocation[1];
    int texLoc2 = faces.at(i).textureLocation[2];
    vec2 tex0 = texturePoints.at(texLoc0 -1);
    vec2 tex1 = texturePoints.at(texLoc1 -1);
    vec2 tex2 = texturePoints.at(texLoc2 -1);


    ModelTriangle newTri = ModelTriangle(pointOne, pointTwo, pointThree, material, tex0, tex1, tex2);
    //cout<<newTri.texture[0].x<<" "<<newTri.texture[0].depth<<" : "<<newTri.texture[1].x<<" "<<newTri.texture[1].y<<" : "<<newTri.texture[2].x<<" "<<newTri.texture[2].y<<endl;
    //cout<<"Point at : "<<newTri.vertices[0].z<<" "<<newTri.vertices[1].z<<" "<<newTri.vertices[2].z<<endl;
    triangles.push_back(newTri);
  }
}


//This will calculate the triangles' location on the canvas
vector<CanvasTriangle>  getCanvas(vector<ModelTriangle> triangles, vec3 cameraLocation, mat3 cameraOrientation, vec2 moveCamera, int textureHeight, int textureWidth){
  vec3 canvasCenter = vec3(0, 0,-40.0f);

  //cameraLocation = cameraLocation * cameraOrientation;

  vector<CanvasTriangle> allCanvasTris;
  vector<CanvasTriangle> allCanvasTrisTEMP;


  vec3 moveIntoFrame = vec3(0,0,0);

  for(int i = 0; i < (signed int)triangles.size(); i++){
    ModelTriangle currentTri = triangles.at(i);  //tri we are about to put into the canvas
    CanvasPoint placementPoints[3];

    for (int p = 0; p < 3; p++){ //loop thought the 3 points in a tri
      vec3 workingPoint = currentTri.vertices[p];
      //Calculate coordinates for the point on the canvas
      float x = (((cameraLocation.x - workingPoint.x) / (cameraLocation.z - canvasCenter.z) ) * -(cameraLocation.z - workingPoint.z ));
      float y = (((cameraLocation.y - workingPoint.y) / (cameraLocation.z - canvasCenter.z) ) * -(cameraLocation.z - workingPoint.z ));
      //The z coordinate, or rather the distance from the camera to the vertice should be calculated, which was attempted below in the commented out line,
      //but we couldn't get the z buffer to work, so we used workingPoint.z because it made it easier to test the texture
      //float z = sqrt(pow((cameraLocation.x - workingPoint.x),2) + pow((cameraLocation.y - workingPoint.y),2) + pow((cameraLocation.z - workingPoint.z),2));


      placementPoints[p] = CanvasPoint(x, y, workingPoint.z);
      //Adjust the coordinates depending on the camera orientation
      vec3 adjustedVector = vec3(x, y, workingPoint.z) * cameraOrientation;
      placementPoints[p] = CanvasPoint(adjustedVector.x, adjustedVector.y, adjustedVector.z);
      placementPoints[p].texturePoint.x = currentTri.texture[p].x * textureWidth;
      placementPoints[p].texturePoint.y = currentTri.texture[p].y * textureHeight;

      if (x < moveIntoFrame.x){
        moveIntoFrame.x = x;
      }

      if (y > moveIntoFrame.y){
        moveIntoFrame.y = y;
      }
    }

    //cout<<placementPoints[0].texturePoint.x<<" "<<placementPoints[0].texturePoint.y<<" : "<<placementPoints[1].texturePoint.x<<" "<<placementPoints[1].texturePoint.y<<" : "<<placementPoints[2].texturePoint.x<<" "<<placementPoints[2].texturePoint.y<<endl;

    allCanvasTrisTEMP.push_back(CanvasTriangle(placementPoints[0], placementPoints[1] ,placementPoints[2], currentTri.colour));

  }

  for(int i = 0; i < (signed int)allCanvasTrisTEMP.size(); i++){
    CanvasTriangle hold = allCanvasTrisTEMP.at(i);
    CanvasPoint p1 = hold.vertices[0];
    CanvasPoint p2 = hold.vertices[1];
    CanvasPoint p3 = hold.vertices[2];

    p1.x = p1.x - moveIntoFrame.x + moveCamera.x;
    p1.y = abs(p1.y - moveIntoFrame.y + moveCamera.y);
    p2.x = p2.x - moveIntoFrame.x + moveCamera.x;
    p2.y = abs(p2.y - moveIntoFrame.y + moveCamera.y);
    p3.x = p3.x - moveIntoFrame.x + moveCamera.x;
    p3.y = abs(p3.y - moveIntoFrame.y + moveCamera.y);

  //  hold = allCanvasTrisTEMP.at(i);
    allCanvasTris.push_back(CanvasTriangle(p1,p2,p3,hold.colour));
    //cout<<p1.texturePoint.x<<" "<<p1.texturePoint.y<<" : "<<p2.texturePoint.x<<" "<<p2.texturePoint.y<<" : "<<p3.texturePoint.x<<" "<<p3.texturePoint.y<<endl;
  }
  return allCanvasTris;
}

void textureLine(CanvasPoint from, CanvasPoint to, vector<Colour> textureMap, int textureHeight, int textureWidth, float buffer[WIDTH][HEIGHT])
{
  //window.clearPixels();
  //cout<<from.texturePoint.x<< " "<<from.texturePoint.y<<endl;
  //cout<<to.texturePoint.x<< " "<<to.texturePoint.y<<endl;
  float xDiff = to.x - from.x;
  float yDiff = to.y - from.y;
  float zDiff = to.depth - from.depth;
  float xTexturePointDiff = to.texturePoint.x - from.texturePoint.x;
  float yTexturePointDiff = to.texturePoint.y - from.texturePoint.y;
  float numberOfSteps1 = std::max(abs(xDiff), abs(yDiff));
  float numberOfSteps2 = std::max(abs(xTexturePointDiff), abs(yTexturePointDiff));
  float numberOfSteps = std::max(abs(numberOfSteps1), abs(numberOfSteps2));
  float xStepSize = xDiff/numberOfSteps;
  float yStepSize = yDiff/numberOfSteps;
  float zStepSize = zDiff/numberOfSteps;
  float xTextureStepSize = xTexturePointDiff/numberOfSteps;
  float yTextureStepSize = yTexturePointDiff/numberOfSteps;
  //uint32_t colourPixel = (255<<24) + (int(colour.red)<<16) + (int( colour.green)<<8) + int(colour.blue);

  //Interpolate between two points to calculate each point's coordinates and the coordinates on the texture map
  for(float i = 0.0; i<numberOfSteps; i++){
    float x = from.x + (xStepSize * i);
    float y = from.y + (yStepSize * i);
    float z = from.depth + (zStepSize * i);
    float xTexture = (round(from.texturePoint.x + (xTextureStepSize * i)));
    float yTexture = (round(from.texturePoint.y + (yTextureStepSize * i)));
    //cout<<xTexture<< " "<<yTexture<<endl;

    //Get colour from the texture map
    Colour textureColour = textureMap.at(xTexture + yTexture *textureWidth);
  //  cout<<textureColour.red<<" "<<textureColour.green<<" "<<textureColour.blue<<endl;
    uint32_t colourPixel = (255<<24) + (int(textureColour.red)<<16) + (int( textureColour.green)<<8) + int( textureColour.blue);
    if(x>=0 && y>=0 && x<=WIDTH && y<=HEIGHT){
      if(buffer[(int)round(x)][(int)round(y)] >= 1/z) {
        buffer[(int)round(x)][(int)round(y)] = 1/z;
        window.setPixelColour(round(x), round(y), colourPixel);
      }
    }
  }
  //cout<<endl;
}

void textureTopTriangle(CanvasTriangle triangle, vector<Colour> textureMap, int textureHeight, int textureWidth, float buffer[WIDTH][HEIGHT])
{
  float Step = 0;
  float PointDiff = triangle.vertices[1].y - triangle.vertices[0].y;
  float TexturePointDiff = triangle.vertices[1].texturePoint.y - triangle.vertices[0].texturePoint.y;
  float numberOfSteps = std::max(abs(PointDiff), abs(TexturePointDiff));
  float StepSize = TexturePointDiff/numberOfSteps;
  //cout<<"Top triangle: "<< TexturePointDiff<<" "<<numberOfSteps<<" "<<StepSize<<endl;

  for(int y = triangle.vertices[0].y; y<=triangle.vertices[1].y; y++){
    //interpolate top to bottom left
    float xleft = triangle.vertices[0].x + (triangle.vertices[1].x - triangle.vertices[0].x)*(y - triangle.vertices[0].y)/(triangle.vertices[1].y - triangle.vertices[0].y);
    //interpolate top to bottom right
    float xright = triangle.vertices[0].x + (triangle.vertices[2].x - triangle.vertices[0].x)*(y - triangle.vertices[0].y)/(triangle.vertices[2].y - triangle.vertices[0].y);
    //do the same for the depth
    float zleft = triangle.vertices[0].depth + (triangle.vertices[1].depth - triangle.vertices[0].depth)*(xleft - triangle.vertices[0].x)/(triangle.vertices[1].x - triangle.vertices[0].x);
    float zright = triangle.vertices[0].depth + (triangle.vertices[2].depth - triangle.vertices[0].depth)*(xright - triangle.vertices[0].x)/(triangle.vertices[2].x - triangle.vertices[0].x);
    CanvasPoint v1 = CanvasPoint(xleft, y, zleft);
    CanvasPoint v2 = CanvasPoint(xright, y, zright);

    //do the same for textures
    float yTexture = triangle.vertices[0].texturePoint.y + Step;
    float textureLeft;
    float textureRight;

    //if texture triangle is the shape of a bottom triangle use the respective interpolation
    if(triangle.vertices[0].texturePoint.y == triangle.vertices[1].texturePoint.y){
      textureLeft = triangle.vertices[0].texturePoint.x + (triangle.vertices[2].texturePoint.x - triangle.vertices[0].texturePoint.x)*(yTexture - triangle.vertices[0].texturePoint.y)/(triangle.vertices[2].texturePoint.y - triangle.vertices[0].texturePoint.y);
      textureRight = triangle.vertices[1].texturePoint.x + (triangle.vertices[2].texturePoint.x - triangle.vertices[1].texturePoint.x)*(yTexture - triangle.vertices[1].texturePoint.y)/(triangle.vertices[2].texturePoint.y - triangle.vertices[1].texturePoint.y);
    }
    //otherwise use the top triangle interpolation
    else {
      textureLeft = triangle.vertices[0].texturePoint.x + (triangle.vertices[1].texturePoint.x - triangle.vertices[0].texturePoint.x)*(yTexture - triangle.vertices[0].texturePoint.y)/(triangle.vertices[1].texturePoint.y - triangle.vertices[0].texturePoint.y);
      textureRight = triangle.vertices[0].texturePoint.x + (triangle.vertices[2].texturePoint.x - triangle.vertices[0].texturePoint.x)*(yTexture - triangle.vertices[0].texturePoint.y)/(triangle.vertices[2].texturePoint.y - triangle.vertices[0].texturePoint.y);

    }
    //draw line by interpolating from left to right
    v1.texturePoint = TexturePoint(textureLeft, yTexture);
    v2.texturePoint = TexturePoint(textureRight, yTexture);
    //cout<<textureLeft<<" "<<textureRight<<" "<<yTexture<<endl;
    //cout<<endl;
    textureLine(v1, v2, textureMap, textureHeight, textureWidth, buffer);
    Step = Step + StepSize;
  }
}

void textureBottomTriangle(CanvasTriangle triangle, vector<Colour> textureMap, int textureHeight, int textureWidth, float buffer[WIDTH][HEIGHT])
{
  float Step = 0;
  float PointDiff = triangle.vertices[2].y - triangle.vertices[0].y;
  float TexturePointDiff = triangle.vertices[2].texturePoint.y - triangle.vertices[0].texturePoint.y;
  float numberOfSteps = std::max(abs(PointDiff), abs(TexturePointDiff));
  float StepSize = TexturePointDiff/numberOfSteps;
  //cout<<triangle.vertices[0].y;
  //cout<<"Bottom triangle: "<<TexturePointDiff<<" "<<numberOfSteps<<" "<<StepSize<<endl;
  //cout<<triangle.vertices[0].texturePoint.x<<" "<<triangle.vertices[1].texturePoint.x<<" "<<triangle.vertices[2].texturePoint.x<<endl;
  for(int y = triangle.vertices[0].y; y<=triangle.vertices[2].y; y++){
    //Interpolate from left to bottom
    float xleft = triangle.vertices[0].x + (triangle.vertices[2].x - triangle.vertices[0].x)*(y - triangle.vertices[0].y)/(triangle.vertices[2].y - triangle.vertices[0].y);
    //Interpolate from right to bottom
    float xright = triangle.vertices[1].x + (triangle.vertices[2].x - triangle.vertices[1].x)*(y - triangle.vertices[1].y)/(triangle.vertices[2].y - triangle.vertices[1].y);
    //do the same for depth
    float zleft = triangle.vertices[0].depth + (triangle.vertices[2].depth - triangle.vertices[0].depth)*(xleft - triangle.vertices[0].x)/(triangle.vertices[2].x - triangle.vertices[0].x);
    float zright = triangle.vertices[1].depth + (triangle.vertices[2].depth - triangle.vertices[1].depth)*(xright - triangle.vertices[1].x)/(triangle.vertices[2].x - triangle.vertices[1].x);
    CanvasPoint v1 = CanvasPoint(xleft, y, zleft);
    CanvasPoint v2 = CanvasPoint(xright, y, zright);

    //do the same for texture
    float yTexture = triangle.vertices[0].texturePoint.y + Step;
    float textureLeft;
    float textureRight;

    //if texture triangle has the shape of a top triangle use the respective interpolation equations
    if(triangle.vertices[1].texturePoint.y == triangle.vertices[2].texturePoint.y){
      textureLeft = triangle.vertices[0].texturePoint.x + (triangle.vertices[1].texturePoint.x - triangle.vertices[0].texturePoint.x)*(yTexture - triangle.vertices[0].texturePoint.y)/(triangle.vertices[1].texturePoint.y - triangle.vertices[0].texturePoint.y);
      textureRight = triangle.vertices[0].texturePoint.x + (triangle.vertices[2].texturePoint.x - triangle.vertices[0].texturePoint.x)*(yTexture - triangle.vertices[0].texturePoint.y)/(triangle.vertices[2].texturePoint.y - triangle.vertices[0].texturePoint.y);
    }
    //otherwise use the bottom triangle interpolation
    else{
      textureLeft = triangle.vertices[0].texturePoint.x + (triangle.vertices[2].texturePoint.x - triangle.vertices[0].texturePoint.x)*(yTexture - triangle.vertices[0].texturePoint.y)/(triangle.vertices[2].texturePoint.y - triangle.vertices[0].texturePoint.y);
      textureRight = triangle.vertices[1].texturePoint.x + (triangle.vertices[2].texturePoint.x - triangle.vertices[1].texturePoint.x)*(yTexture - triangle.vertices[1].texturePoint.y)/(triangle.vertices[2].texturePoint.y - triangle.vertices[1].texturePoint.y);
    }

    //interpolate from left to right
    v1.texturePoint = TexturePoint(textureLeft, yTexture);
    v2.texturePoint = TexturePoint(textureRight, yTexture);
    //cout<<textureLeft<<" "<<textureRight<<" "<<yTexture<<endl;

    textureLine(v1, v2, textureMap, textureHeight, textureWidth, buffer);
    Step = Step + StepSize;
  }
}

void textureTriangle(CanvasTriangle triangle,  vector<Colour> textureMap, int textureHeight, int textureWidth, float buffer[WIDTH][HEIGHT])
{

  triangle = sortVertices(triangle);     //sort triangle points

  //cout<<"Triangle z coordinates: "<<triangle.vertices[0].depth<<" "<<triangle.vertices[1].depth<<" "<<triangle.vertices[2].depth<<endl;
  //If the triangle is in the shape of a top triangle we don't need to split it
  if(triangle.vertices[1].y == triangle.vertices[2].y){
    if(triangle.vertices[1].x > triangle.vertices[2].x){
      CanvasPoint hold = triangle.vertices[1];
      triangle.vertices[1] = triangle.vertices[2];
      triangle.vertices[2] = hold;
    }
    //store texture points into a CanvasTriangle
    CanvasPoint tex0 = CanvasPoint(triangle.vertices[0].texturePoint.x, triangle.vertices[0].texturePoint.y);
    CanvasPoint tex1 = CanvasPoint(triangle.vertices[1].texturePoint.x, triangle.vertices[1].texturePoint.y);
    CanvasPoint tex2 = CanvasPoint(triangle.vertices[2].texturePoint.x, triangle.vertices[2].texturePoint.y);
    CanvasTriangle textureTriangle = CanvasTriangle(tex0, tex1, tex2);
    textureTriangle = sortVertices(textureTriangle);        //sort texture points

    //Set the corresponding texture points after sorting them
    triangle.vertices[0].texturePoint = TexturePoint(textureTriangle.vertices[0].x, textureTriangle.vertices[0].y);
    triangle.vertices[1].texturePoint = TexturePoint(textureTriangle.vertices[1].x, textureTriangle.vertices[1].y);
    triangle.vertices[2].texturePoint = TexturePoint(textureTriangle.vertices[2].x, textureTriangle.vertices[2].y);

    //cout<<"Top triangle: "<<triangle.vertices[0].texturePoint.y<<" "<<triangle.vertices[1].texturePoint.y<<" "<<triangle.vertices[2].texturePoint.y<<endl;
    if(triangle.vertices[0].texturePoint.y != triangle.vertices[1].texturePoint.y || triangle.vertices[0].texturePoint.y != triangle.vertices[2].texturePoint.y || triangle.vertices[1].texturePoint.y != triangle.vertices[2].texturePoint.y)
      textureTopTriangle(triangle, textureMap, textureHeight, textureWidth, buffer);
  }
  //If the triangle is in the shape of a bottom triangle we also don't need to split it
  else if(triangle.vertices[0].y == triangle.vertices[1].y){
    if(triangle.vertices[0].x > triangle.vertices[1].x){
      CanvasPoint hold = triangle.vertices[1];
      triangle.vertices[1] = triangle.vertices[0];
      triangle.vertices[0] = hold;
    }

    //store texture points into a CanvasTriangle
    CanvasPoint tex0 = CanvasPoint(triangle.vertices[0].texturePoint.x, triangle.vertices[0].texturePoint.y);
    CanvasPoint tex1 = CanvasPoint(triangle.vertices[1].texturePoint.x, triangle.vertices[1].texturePoint.y);
    CanvasPoint tex2 = CanvasPoint(triangle.vertices[2].texturePoint.x, triangle.vertices[2].texturePoint.y);
    CanvasTriangle textureTriangle = CanvasTriangle(tex0, tex1, tex2);
    textureTriangle = sortVertices(textureTriangle);        //sort texture points

    //Set the corresponding texture points after sorting them
    triangle.vertices[0].texturePoint = TexturePoint(textureTriangle.vertices[0].x, textureTriangle.vertices[0].y);
    triangle.vertices[1].texturePoint = TexturePoint(textureTriangle.vertices[1].x, textureTriangle.vertices[1].y);
    triangle.vertices[2].texturePoint = TexturePoint(textureTriangle.vertices[2].x, textureTriangle.vertices[2].y);

    //cout<<"Bottom triangle: "<<triangle.vertices[0].texturePoint.y<<" "<<triangle.vertices[1].texturePoint.y<<" "<<triangle.vertices[2].texturePoint.y<<endl;
    if(triangle.vertices[0].texturePoint.y != triangle.vertices[1].texturePoint.y || triangle.vertices[0].texturePoint.y != triangle.vertices[2].texturePoint.y || triangle.vertices[1].texturePoint.y != triangle.vertices[2].texturePoint.y)
      textureBottomTriangle(triangle, textureMap, textureHeight, textureWidth, buffer);
  }
  else{
    //store texture points into a CanvasTriangle
    CanvasPoint tex0 = CanvasPoint(triangle.vertices[0].texturePoint.x, triangle.vertices[0].texturePoint.y);
    CanvasPoint tex1 = CanvasPoint(triangle.vertices[1].texturePoint.x, triangle.vertices[1].texturePoint.y);
    CanvasPoint tex2 = CanvasPoint(triangle.vertices[2].texturePoint.x, triangle.vertices[2].texturePoint.y);
    CanvasTriangle textureTriangle = CanvasTriangle(tex0, tex1, tex2);
    textureTriangle = sortVertices(textureTriangle);        //sort texture points

    //Set the corresponding texture points after sorting them
    triangle.vertices[0].texturePoint = TexturePoint(textureTriangle.vertices[0].x, textureTriangle.vertices[0].y);
    triangle.vertices[1].texturePoint = TexturePoint(textureTriangle.vertices[1].x, textureTriangle.vertices[1].y);
    triangle.vertices[2].texturePoint = TexturePoint(textureTriangle.vertices[2].x, textureTriangle.vertices[2].y);

    //Create a new vertice
    CanvasPoint vertice = createVertice(triangle.vertices[0],triangle.vertices[1],triangle.vertices[2]);

    //Use the new vertice to split triangle into two
    CanvasTriangle triangle1 = CanvasTriangle(triangle.vertices[0],vertice,triangle.vertices[1]);
    //put the vertices in a specific order
    if(triangle1.vertices[1].x > triangle1.vertices[2].x){
      CanvasPoint hold = triangle1.vertices[1];
      triangle1.vertices[1] = triangle1.vertices[2];
      triangle1.vertices[2] = hold;

      CanvasPoint tex0 = CanvasPoint(triangle1.vertices[0].texturePoint.x, triangle1.vertices[0].texturePoint.y);
      CanvasPoint tex1 = CanvasPoint(triangle1.vertices[1].texturePoint.x, triangle1.vertices[1].texturePoint.y);
      CanvasPoint tex2 = CanvasPoint(triangle1.vertices[2].texturePoint.x, triangle1.vertices[2].texturePoint.y);
      CanvasTriangle textureTriangle = CanvasTriangle(tex0, tex1, tex2);
      textureTriangle = sortVertices(textureTriangle);        //sort texture points

      //Set the corresponding texture points after sorting them
      triangle1.vertices[0].texturePoint = TexturePoint(textureTriangle.vertices[0].x, textureTriangle.vertices[0].y);
      triangle1.vertices[1].texturePoint = TexturePoint(textureTriangle.vertices[1].x, textureTriangle.vertices[1].y);
      triangle1.vertices[2].texturePoint = TexturePoint(textureTriangle.vertices[2].x, textureTriangle.vertices[2].y);

    }
    CanvasTriangle triangle2 = CanvasTriangle(triangle.vertices[1],vertice,triangle.vertices[2]);
    //put the vertices in a specific order
    if(triangle2.vertices[0].x > triangle2.vertices[1].x){
      CanvasPoint hold = triangle2.vertices[1];
      triangle2.vertices[1] = triangle2.vertices[0];
      triangle2.vertices[0] = hold;

      CanvasPoint tex0 = CanvasPoint(triangle2.vertices[0].texturePoint.x, triangle2.vertices[0].texturePoint.y);
      CanvasPoint tex1 = CanvasPoint(triangle2.vertices[1].texturePoint.x, triangle2.vertices[1].texturePoint.y);
      CanvasPoint tex2 = CanvasPoint(triangle2.vertices[2].texturePoint.x, triangle2.vertices[2].texturePoint.y);
      CanvasTriangle textureTriangle = CanvasTriangle(tex0, tex1, tex2);
      textureTriangle = sortVertices(textureTriangle);        //sort texture points

      //Set the corresponding texture points after sorting them
      triangle2.vertices[0].texturePoint = TexturePoint(textureTriangle.vertices[0].x, textureTriangle.vertices[0].y);
      triangle2.vertices[1].texturePoint = TexturePoint(textureTriangle.vertices[1].x, textureTriangle.vertices[1].y);
      triangle2.vertices[2].texturePoint = TexturePoint(textureTriangle.vertices[2].x, textureTriangle.vertices[2].y);

    }

    //cout<<"Split Top triangle: "<<triangle1.vertices[0].texturePoint.y<<" "<<triangle1.vertices[1].texturePoint.y<<" "<<triangle1.vertices[2].texturePoint.y<<endl;

    if(triangle1.vertices[0].texturePoint.y != triangle1.vertices[1].texturePoint.y || triangle1.vertices[0].texturePoint.y != triangle1.vertices[2].texturePoint.y || triangle1.vertices[1].texturePoint.y != triangle1.vertices[2].texturePoint.y)
      textureTopTriangle(triangle1, textureMap, textureHeight, textureWidth, buffer);
    //cout<<"Split Bottom triangle: "<<triangle2.vertices[0].texturePoint.y<<" "<<triangle2.vertices[1].texturePoint.y<<" "<<triangle2.vertices[2].texturePoint.y<<endl;
    if(triangle2.vertices[0].texturePoint.y != triangle2.vertices[1].texturePoint.y || triangle2.vertices[0].texturePoint.y != triangle2.vertices[2].texturePoint.y || triangle2.vertices[1].texturePoint.y != triangle2.vertices[2].texturePoint.y)
      textureBottomTriangle(triangle2, textureMap,  textureHeight, textureWidth, buffer);
  }

}


void drawWireframe(vector<CanvasTriangle> allTri, float buffer[WIDTH][HEIGHT]){
  Colour colour = Colour(255,0,0);
  for(int i = 0; i < (signed int)allTri.size(); i++){
    drawTriangle(allTri.at(i), colour, buffer);
  }
}

void rasterizer(vector<CanvasTriangle> allTri, float buffer[WIDTH][HEIGHT]){
  for(int i = 0; i < (signed int)allTri.size(); i++){
    drawFilledTriangle(allTri.at(i), allTri.at(i).colour, buffer);
  }
}

void textureMapping(vector<CanvasTriangle>allTri, vector<Colour>textureMap, int textureHeight, int textureWidth, float buffer[WIDTH][HEIGHT]){
  cout<<"Starting texture"<<endl;
  for(int i = 0; i < (signed int)allTri.size(); i++){
    textureTriangle(allTri.at(i), textureMap, textureHeight, textureWidth, buffer);
  }
}

void savePPMfile(unsigned int height, unsigned int width, unsigned int colour_max, int &ppm){
  ppm += 1;
  string name;
  if(ppm<10) name = "0000" + std::to_string(ppm) + ".ppm";
  else if(ppm<100) name = "000" + std::to_string(ppm) + ".ppm";
  else if(ppm<1000) name = "00" + std::to_string(ppm) + ".ppm";
  else if(ppm<10000) name = "0" + std::to_string(ppm) + ".ppm";
  else if(ppm<100000) name = std::to_string(ppm) + ".ppm";
  const char * filename = name.c_str();
  FILE *f = fopen(filename, "w");
  fprintf(f, "P6\n%d %d\n%d\n", width, height, colour_max);
  for(int y=0; y<(signed int)height; y++){
    for(int x=0; x<(signed int)width; x++){
      uint32_t pixel = window.getPixelColour(x,y);
      Uint8 red = Uint8(pixel>>16);
      Uint8 green = Uint8(pixel>>8);
      Uint8 blue = Uint8(pixel);

      fprintf(f, "%c%c%c", red, green, blue);
    }
  }
  fclose(f);
}

int main(int argc, char* argv[])
{
  vector<Colour> palette;
  vector<Colour> textureMap;
  vector<ModelTriangle> triangles;
  int textureHeight, textureWidth;
  bool filled = false, texture = false;
  int ppm = 0;
  float buffer[WIDTH][HEIGHT];

  for(int i=0; i<WIDTH; i++)
    for(int j=0; j<HEIGHT; j++)
      buffer[i][j] = std::numeric_limits<float>::infinity();

  //readMtlFile(palette);
  readObjFile(triangles, palette, textureMap, textureHeight, textureWidth);
  vec3 cameraLocation = vec3(0.0f, 30.0f, -46.0f);
  vec2 moveCamera = vec2(0.0f, 0.0f);

  mat3 cameraOrientation(vec3(1.0, 0.0, 0.0), vec3(0.0, 1.0, 0.0), vec3(0.0, 0.0, 1.0));
  vector<CanvasTriangle> allCanvasTris =  getCanvas(triangles, cameraLocation, cameraOrientation, moveCamera, textureHeight, textureWidth);


/*  for(int i = 0; i < (signed int)allCanvasTris.size(); i++){
    cout << "Triangle on the canvas: " << allCanvasTris.at(i) << endl;
  }*/

  //rasterizer(allCanvasTris, buffer);
  //textureMapping(allCanvasTris, textureMap, textureHeight, textureWidth, buffer);
  drawWireframe(allCanvasTris, buffer);

  cout<<"Press LEFT/ RIGHT/ UP/ DOWN to move the camera"<<endl;
  cout<<"Press W to zoom in"<<endl;
  cout<<"Press S to zoom out"<<endl;
  cout<<"Press P to save the image as a PPM file"<<endl;
  cout<<"Press X to rotate camera about X axis"<<endl;
  cout<<"Press Y to rotate camera about Y axis"<<endl;
  cout<<"Press Z to rotate camera about Z axis"<<endl;
  cout<<"Press F to fill in the traingles"<<endl;
  cout<<"Press T for texture"<<endl;

  SDL_Event event;
  while(true)
    {// We MUST poll for events - otherwise the window will freeze !
    if(window.pollForInputEvents(&event)) {
      handleEvent(event, cameraLocation, cameraOrientation, moveCamera, filled, texture, ppm);
      window.clearPixels();
      allCanvasTris =  getCanvas(triangles, cameraLocation, cameraOrientation, moveCamera, textureHeight, textureWidth);
      for(int i=0; i<WIDTH; i++)
        for(int j=0; j<HEIGHT; j++)
          buffer[i][j] = std::numeric_limits<float>::infinity();
      if (filled) rasterizer(allCanvasTris, buffer);
      if (texture) textureMapping(allCanvasTris, textureMap, textureHeight, textureWidth, buffer);

      drawWireframe(allCanvasTris, buffer);
    }
    update();
    // Need to render the frame at the end, or nothing actually gets shown on the screen !
    window.renderFrame();}


  return 0;
}
