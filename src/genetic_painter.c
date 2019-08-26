/*
 * genetic_painter.c
 * Copyright (C) 2019 transpalette <theo.morales.fr@gmail.com>
 *
 * Distributed under terms of the MIT license.
 */

#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <GL/glew.h>
#include <GL/glut.h>
#include <GLFW/glfw3.h>

#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image.h"
#include "stb_image_write.h"

#define GENERATION_SIZE 30
#define MUTATION_RATE 0.01
#define MUTATION_AMOUNT 0.1
#define SELECTION_CUTOFF 0.25
#define POLYGONS 10
#define RENDER_W 200
#define RENDER_H 200

#define random_coord() (rand()/(float)RAND_MAX) * 2.0f - 1.0f
#define rand_a_b(a, b) (((float)rand())/(float)RAND_MAX)*(b-a)+a
#define max(a,b) \
  ({ __typeof__ (a) _a = (a); \
   __typeof__ (b) _b = (b); \
   _a > _b ? _a : _b; })
#define min(a,b) \
  ({ __typeof__ (a) _a = (a); \
   __typeof__ (b) _b = (b); \
   _a < _b ? _a : _b; })

typedef struct
{
  float r, g, b, a;
} Color;

typedef struct
{
  float x, y;
} Vertex;

typedef struct
{
  Vertex vertices[3];
  Color color;
} Triangle;

typedef struct
{
  Triangle genes[POLYGONS];
  float fitness;
  long score;
} Chromosome;

typedef struct
{
  long minScore, maxScore;
  long meanScore;
  float meanFitness;
  Chromosome chromosomes[GENERATION_SIZE];
  Chromosome *fittest;
} Generation;

GLuint program;
GLFWwindow *window;
unsigned char *image;
int imageWidth, imageHeight, imageChannels;


void make_color(Color *c, float r, float g, float b, float a)
{
  c->r = r;
  c->g = g;
  c->b = b;
  c->a = a;
}


void make_vertex(Vertex *v, float x, float y)
{
  v->x = x;
  v->y = y;
}


GLuint make_vao(Triangle *triangles, int length)
{
  GLuint vbo, vao, colorBuffer;

  float points[length*3*3]; // 3 coordinates * 3 vertices
  float colors[length*4*3]; // 4 color channels (RGBA) * 3 vertices

  for (int i = 0, n = 0; i < length*3*3; i += 9, n++) {
    for (int j = 0; j < 9; j += 3) {
      points[i+j] = triangles[n].vertices[j/3].x;
      points[i+(j+1)] = triangles[n].vertices[j/3].y;
      points[i+(j+2)] = 0;
    }
  }

  for (int i = 0, n = 0; i < length*3*4; i += 12, n++) {
    for (int j = 0; j < 12; j += 4) {
      colors[i+j] = triangles[n].color.r;
      colors[i+(j+1)] = triangles[n].color.g;
      colors[i+(j+2)] = triangles[n].color.b;
      colors[i+(j+3)] = triangles[n].color.a;
    }
  }

  vbo = vao = colorBuffer = 0;
  glGenBuffers(1, &vbo);
  glBindBuffer(GL_ARRAY_BUFFER, vbo);
  glBufferData(GL_ARRAY_BUFFER, length * 3 * 3 * sizeof(float), points, GL_STATIC_DRAW);

  glGenBuffers(1, &colorBuffer);
  glBindBuffer(GL_ARRAY_BUFFER, colorBuffer);
  glBufferData(GL_ARRAY_BUFFER, length * 4 * 3 * sizeof(float), colors, GL_STATIC_DRAW);

  glGenVertexArrays(1, &vao);
  glBindVertexArray(vao);
  glEnableVertexAttribArray(0);
  glBindBuffer(GL_ARRAY_BUFFER, vbo);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, NULL);

  glEnableVertexAttribArray(1);
  glBindBuffer(GL_ARRAY_BUFFER, colorBuffer);
  glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, 0, NULL);

  return vao;
}


void draw(Chromosome *c)
{
  GLuint vao = make_vao(c->genes, POLYGONS);


  // wipe the drawing surface clear
  glClearColor(0, 0, 0, 0);
  glClear(GL_COLOR_BUFFER_BIT);
  glUseProgram(program);
  glBindVertexArray(vao);
  // draw points 0-3 from the currently bound VAO with current in-use shader
  glDrawArrays(GL_TRIANGLES, 0, 3 * POLYGONS);
  // update other events like input handling
  glfwPollEvents();
  // put the stuff we've been drawing onto the display
  glfwSwapBuffers(window);
}


void save_jpg(Generation *gen, int n)
{
  unsigned char frame[3*RENDER_W*RENDER_H];
  char name[100];
  sprintf(name, "generation-%d_fittest_score-%ld.png", n, gen->fittest->score);
  draw(gen->fittest);
  glReadPixels(0, 0, RENDER_W, RENDER_H, GL_RGB, GL_UNSIGNED_BYTE, &frame[0]);
  stbi_write_jpg(name, RENDER_W, RENDER_H, 3, frame, 100);
  printf("Saved the fittest chromosome as '%s'", name);
}


void compute_score(Chromosome *c)
{
  long score = 0;
  unsigned char frame[3*RENDER_W*RENDER_H];

  draw(c);
  glReadPixels(0, 0, RENDER_W, RENDER_H, GL_RGB, GL_UNSIGNED_BYTE, &frame[0]);
  for (int i = 0; i < imageWidth*imageHeight*3; i+=3) {
	score +=
	  (frame[i]-image[i])*(frame[i]-image[i])
	  +(frame[i+1]-image[i+1])*(frame[i+1]-image[i+1])
	  +(frame[i+2]-image[i+2])*(frame[i+2]-image[i+2]);
  }

  c->score = score;
}


void compute_fitness(Generation *g, Chromosome *c)
{
  if (c->score == g->maxScore) {
	c->fitness = 0.99;
  	return;
  }

  float normalizedScore = (c->score - g->minScore) / (g->maxScore - g->minScore);

  c->fitness =  abs(normalizedScore - 1);
}


void make_chromosome(Chromosome *c)
{
  for (int i = 0; i < POLYGONS; i++) {
  	for (int j = 0; j < 3; j++) {
	  make_vertex(&c->genes[i].vertices[j], random_coord(), random_coord());
	}
	make_color(&c->genes[i].color, rand_a_b(0, 1), rand_a_b(0, 1), rand_a_b(0, 1), .5f);
  }
}


Chromosome crossover(Chromosome *candidates, int nbCandidates, int index)
{
  Chromosome kid;
  int firstPick, secondPick, crossoverPoint;
  short picked = 0;

  for (int i = 0; i < nbCandidates; i++) {
  	float wheel = rand() / (float) RAND_MAX;
  	if (wheel <= candidates[i].fitness) {
	  firstPick = i;
	  picked = 1;
	  break;
	}
  }

  if (!picked) {
  	firstPick = rand_a_b(0, nbCandidates);
  }

  picked = 0;
  for (int i = 0; i < nbCandidates; i++) {
  	float wheel = rand() / (float) RAND_MAX;
  	if (i != firstPick && wheel <= candidates[i].fitness) {
	  secondPick = i;
	  picked = 1;
	  break;
	}
  }

  if (!picked) {
  	secondPick = rand_a_b(0, nbCandidates);
  }

  crossoverPoint = rand_a_b(0, POLYGONS);
  for (int i = 0; i < POLYGONS; i++) {
	kid.genes[i] = i <= crossoverPoint ? candidates[firstPick].genes[i]
	  : candidates[secondPick].genes[i];
  }

  return kid;
}

void mutate(Chromosome *c)
{
  float wheel;
  for (int i = 0; i < POLYGONS; i++) {
	for (int v = 0; v < 3; v++) {
	  wheel = rand() / (float) RAND_MAX;
	  if (wheel <= MUTATION_RATE) {
		c->genes[i].vertices[v].x += random_coord()*MUTATION_AMOUNT*2 - MUTATION_AMOUNT;
		if (c->genes[i].vertices[v].x < 0) {
		  c->genes[i].vertices[v].x = max(c->genes[i].vertices[v].x, -1);
		} else {
		  c->genes[i].vertices[v].x = min(c->genes[i].vertices[v].x, 1);
		}
	  }
	  wheel = rand() / (float) RAND_MAX;
	  if (wheel <= MUTATION_RATE) {
		c->genes[i].vertices[v].y += random_coord()*MUTATION_AMOUNT*2 - MUTATION_AMOUNT;
		if (c->genes[i].vertices[v].y < 0) {
		  c->genes[i].vertices[v].y = max(c->genes[i].vertices[v].y, -1);
		} else {
		  c->genes[i].vertices[v].y = min(c->genes[i].vertices[v].y, 1);
		}
	  }
	}
	for (int j = 0; j < 3; j++) {
	  wheel = rand() / (float) RAND_MAX;
	  if (wheel <= MUTATION_RATE) {
	  	switch(j) {
		  case 0:
			c->genes[i].color.r += rand_a_b(0, 1)*MUTATION_AMOUNT*2 - MUTATION_AMOUNT;
			if (c->genes[i].color.r < 0) {
			  c->genes[i].color.r = max(c->genes[i].color.r, 0);
			} else {
			  c->genes[i].color.r = min(c->genes[i].color.r, 1);
			}
			break;
		  case 1:
			c->genes[i].color.g += rand_a_b(0, 1)*MUTATION_AMOUNT*2 - MUTATION_AMOUNT;
			if (c->genes[i].color.g < 0) {
			  c->genes[i].color.g = max(c->genes[i].color.g, 0);
			} else {
			  c->genes[i].color.g = min(c->genes[i].color.g, 1);
			}
		  	break;
		  case 2:
			c->genes[i].color.b += rand_a_b(0, 1)*MUTATION_AMOUNT*2 - MUTATION_AMOUNT;
			if (c->genes[i].color.b < 0) {
			  c->genes[i].color.b = max(c->genes[i].color.b, 0);
			} else {
			  c->genes[i].color.b = min(c->genes[i].color.b, 1);
			}
		  case 3:
			c->genes[i].color.a += rand_a_b(0, 1)*MUTATION_AMOUNT*2 - MUTATION_AMOUNT;
			if (c->genes[i].color.a < 0) {
			  c->genes[i].color.a = max(c->genes[i].color.a, 0);
			} else {
			  c->genes[i].color.a = min(c->genes[i].color.a, 1);
			  }
		}
	  }
	}
  }
}

GLFWwindow* init_opengl()
{
  if (!glfwInit()) {
    fprintf(stderr, "ERROR: could not start GLFW3\n");
	exit(1);
  }

  glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

  GLFWwindow* window = glfwCreateWindow(
  	  RENDER_W, RENDER_H, "MonaLisa Genetic Solver", NULL, NULL);

  if (!window) {
    fprintf(stderr, "ERROR: could not open window with GLFW3\n");
    glfwTerminate();
    exit(1);
  }

  glfwMakeContextCurrent(window);
  // start GLEW extension handler
  glewExperimental = GL_TRUE;
  glewInit();

  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glEnable(GL_BLEND);

  return window;
}


GLuint build_program()
{
  const char* vertex_shader =
	"#version 410\n"
	"layout (location = 0) in vec3 vp;"
	"layout (location = 1) in vec4 vertexColor;"
	"out vec4 fragmentColor;"
	"void main() {"
	"  gl_Position = vec4(vp, 1.0);"
	"  fragmentColor = vertexColor;"
	"}";

  const char* fragment_shader =
	"#version 410\n"
	"out vec4 frag_color;"
	"in vec4 fragmentColor;"
	"void main() {"
	"  frag_color = fragmentColor;"
	"}";

  GLuint vs = glCreateShader(GL_VERTEX_SHADER);
  glShaderSource(vs, 1, &vertex_shader, NULL);
  glCompileShader(vs);
  GLuint fs = glCreateShader(GL_FRAGMENT_SHADER);
  glShaderSource(fs, 1, &fragment_shader, NULL);
  glCompileShader(fs);

  GLuint shaderProgram = glCreateProgram();
  glAttachShader(shaderProgram, fs);
  glAttachShader(shaderProgram, vs);
  glLinkProgram(shaderProgram);


  return shaderProgram;
}


void compute_metrics(Generation *gen)
{
  long meanScore = 0;
  int fittest = 0;
  float meanFitness = 0;
  gen->minScore = 0;
  gen->maxScore = 0;

  for (int i = 0; i < GENERATION_SIZE; i++) {
  	compute_score(&gen->chromosomes[i]);
	if (gen->chromosomes[i].score > gen->maxScore) {
	  gen->maxScore = gen->chromosomes[i].score;
	}
	if (gen->minScore == 0 || gen->chromosomes[i].score < gen->minScore) {
	  gen->minScore = gen->chromosomes[i].score;
	}

	meanScore += gen->chromosomes[i].score;
  }

  /* Compute each chromosome's fitness now that maxScore & minScore are
   * computed */
  for (int i = 0; i < GENERATION_SIZE; i++) {
  	compute_fitness(gen, &gen->chromosomes[i]);
  	if (gen->chromosomes[i].fitness > fittest) {
	  fittest = i;
	}
  	meanFitness += gen->chromosomes[i].fitness;
  }

  gen->meanScore = meanScore / GENERATION_SIZE;
  gen->meanFitness = meanFitness / GENERATION_SIZE;
  gen->fittest = &gen->chromosomes[fittest];
}


int comp_chromosomes(const void *el1, const void *el2)
{
  Chromosome f = *(Chromosome*)el1;
  Chromosome s = *(Chromosome*)el2;

  return f.score - s.score;
}


/* arrange the N elements of ARRAY in random order.
 * Only effective if N is much smaller than RAND_MAX;
 * if this may not be the case, use a better random
 * number generator. */
static void shuffle(void *array, size_t n, size_t size) {
    char tmp[size];
    char *arr = array;
    size_t stride = size * sizeof(char);

    if (n > 1) {
        size_t i;
        for (i = 0; i < n - 1; ++i) {
            size_t rnd = (size_t) rand();
            size_t j = i + rnd / (RAND_MAX / (n - i) + 1);

            memcpy(tmp, arr + j * stride, size);
            memcpy(arr + j * stride, arr + i * stride, size);
            memcpy(arr + i * stride, tmp, size);
        }
    }
}


/* Selects the generation's chromosomes from the given population/generation,
 * by randomly picking them, but with a higher probability for the fittest
 * chromosomes
 */
void select_from_generation(Generation *gen, int amount, Chromosome *selection_dest)
{
  short picked, skip;
  int selected[amount];

  memset(selected, -1, amount-1);
  for (int i = 0; i < amount; i++) {
	picked = 0;
	for (int j = 0; j < GENERATION_SIZE; j++) {
	  float wheel = rand_a_b(0, 1);
	  if (wheel <= gen->chromosomes[j].fitness) {
	  	skip = 0;
	  	for (int k = 0; k < i; k++) {
		  if (selected[k] == j) {
		  	skip = 1;
		  	break;
		  }
		}
		if (!skip) {
		  selection_dest[i] = gen->chromosomes[j];
		  mutate(&selection_dest[i]);
		  selected[i] = j;
		  picked = 1;
		  break;
		}
	  }
	}

	if (!picked) {
	  selection_dest[i] = gen->chromosomes[(int)rand_a_b(0, GENERATION_SIZE)];
	  mutate(&selection_dest[i]);
	}
  }
}


void evolve(Generation *gen, int index)
{
  int nbOffsprings = (int)(SELECTION_CUTOFF*(float)GENERATION_SIZE);
  Chromosome offsprings[nbOffsprings];
  Chromosome selection[GENERATION_SIZE-nbOffsprings];

  compute_metrics(gen);
  printf("[*] Generation %d:\n\t-> Average score: %ld\n\t-> Fittest's score: %ld\n\n",
	  index, gen->meanScore, gen->fittest->score);

  qsort(gen->chromosomes, GENERATION_SIZE, sizeof(Chromosome), comp_chromosomes);

  for (int i = 0; i < nbOffsprings; i++) {
	offsprings[i] = crossover(gen->chromosomes, GENERATION_SIZE, i);
  }

  shuffle(gen->chromosomes, GENERATION_SIZE, sizeof(Chromosome));
  select_from_generation(gen, GENERATION_SIZE-nbOffsprings, selection);
  for (int i = 0; i < GENERATION_SIZE; i++) {
	gen->chromosomes[i] = i < nbOffsprings ? offsprings[i] : selection[i];
  }
}


Generation *make_random_generation()
{
  Generation *gen = (Generation*) malloc(sizeof(Generation));
  for (int n = 0; n < GENERATION_SIZE; n++)
  	make_chromosome(&gen->chromosomes[n]);

  return gen;
}

int main(int argc, char **argv) {
  if (argc < 2) {
  	printf("Usage: %s <picture.png>\n", argv[0]);
	return 1;
  }

  srand((unsigned) time(NULL));
  window = init_opengl();
  program = build_program();

  image = stbi_load(argv[1], &imageWidth, &imageHeight, &imageChannels, 0);
  Generation *gen = make_random_generation();
  int n = 0;

  while(!glfwWindowShouldClose(window)) {
	evolve(gen, n++);
	if (n % 100 == 0) {
	  save_jpg(gen, n);
	}
  }

  save_jpg(gen, n);
  // close GL context and any other GLFW resources
  glfwTerminate();
  stbi_image_free(image);

  return 0;
}
