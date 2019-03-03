#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <pthread.h>
#include <stdint.h>

#define POS_X  0
#define POS_Y 1
#define MASS 2

const double epsilon_0 = 0.001;

//declare node structure (X)
//maybe good idea is to compact it! (ie no padding)
typedef struct tree_node{
  int depth;
  int body_id;
  int is_used;
  double x_lim; //x division (middle point in x)
  double y_lim; // y division (middle point in y)
  double width; //width of the box
  double cm_x; // = 0, this points to 0 and
  double cm_y; //center of mas s of the quadrant
  double tot_mass; // = 0, total mass of the quadrant
  struct tree_node *left_down; //Q3 child
  struct tree_node *left_up; //Q2 child
  struct tree_node *right_down; //Q4 child
  struct tree_node *right_up; //Q1 child
}node_t;

void print_qtree(node_t* node){
  if(node == NULL){printf("Tree is empty \n"); return;}

  if(node != NULL){
    printf("Depth: %d, id:%d, total mass: %lf, cm: (%lf, %lf) limits: (%lf,%lf)\t",
    node->depth,
    node->body_id,
    node->tot_mass,
    node->cm_x, node->cm_y,
    node->x_lim, node->y_lim);
  }

  if (node->right_up != NULL) printf("RU : total mass %lf,", node->right_up->tot_mass);
  if (node->left_up != NULL) printf("LU: total mass: %lf,", node->left_up->tot_mass);
  if (node->left_down != NULL) printf("LD: total mass: %lf,", node->left_down->tot_mass);
  if (node->right_down != NULL) printf("RD: total mass; %lf,", node->right_down->tot_mass);
  printf("\n");

  if (node->right_up != NULL) print_qtree(node->right_up);
  if (node->left_up != NULL) print_qtree(node->left_up);
  if (node->left_down != NULL) print_qtree(node->left_down);
  if (node->right_down != NULL) print_qtree(node->right_down);
  return;
}
void create_children(node_t* node, double * pow_2){

  if(node->right_up != NULL){printf("L178, create_children. Something went wrong RU children already exist"); return;}
  //create right up
  //printf("Creating RU\n");
  node->right_up=(node_t*)malloc(sizeof(node_t));
  node->right_up->is_used = 0;
  node->right_up->depth = node->depth +1;
  node->right_up->body_id = -1;
  node->right_up->x_lim = node->x_lim + pow_2[node->right_up->depth];
  node->right_up->y_lim =  node->y_lim + pow_2[node->right_up->depth];
  node->right_up-> width = pow_2[node->depth];
  node->right_up-> cm_x = 0;
  node->right_up-> cm_y = 0;
  node->right_up-> tot_mass = -1.0;
  node->right_up->left_down = NULL;
  node->right_up->left_up = NULL;
  node->right_up->right_down = NULL;
  node->right_up->right_up = NULL;

  if(node->right_down != NULL){printf("L178, create_children. Something went wrong RD children already exist"); return;}
  //create right down
  //printf("Creating RD\n");
  node->right_down=(node_t*)malloc(sizeof(node_t));
  node->right_down->is_used=0;
  node->right_down->depth = node->depth +1;
  node->right_down->body_id = -1;
  node->right_down->x_lim = node->x_lim + pow_2[node->right_down->depth];
  node->right_down->y_lim =  node->y_lim - pow_2[node->right_down->depth];
  node->right_down-> width = pow_2[node->depth];
  node->right_down-> cm_x = 0;
  node->right_down-> cm_y = 0;
  node->right_down-> tot_mass = -1.0;
  node->right_down->left_down = NULL;
  node->right_down->left_up = NULL;
  node->right_down->right_down = NULL;
  node->right_down->right_up = NULL;


  if(node->left_up != NULL){printf("L178, create_children. Something went wrong RU children already exist"); return;}
  //create left_up
  //printf("Creating LU\n");
  node->left_up=(node_t*)malloc(sizeof(node_t));
  node->left_up->is_used = 0;
  node->left_up->depth = node->depth +1;
  node->left_up->body_id = -1;
  node->left_up->x_lim = node->x_lim - pow_2[node->left_up->depth];
  node->left_up->y_lim =  node->y_lim + pow_2[node->left_up->depth];
  node->left_up-> width = pow_2[node->depth];
  node->left_up-> cm_x = 0;
  node->left_up-> cm_y = 0;
  node->left_up-> tot_mass = -1.0;
  node->left_up->left_down = NULL;
  node->left_up->left_up = NULL;
  node->left_up->right_down = NULL;
  node->left_up->right_up = NULL;


  if(node->left_down != NULL){printf("L178, create_children. Something went wrong RU children already exist"); return;}
  //create left down
  //printf("Creating LD\n");
  node->left_down=(node_t*)malloc(sizeof(node_t));
  node->left_down->is_used=0;
  node->left_down->depth = node->depth +1;
  node->left_down->body_id = -1;
  node->left_down->x_lim = node->x_lim - pow_2[node->left_down->depth];
  node->left_down->y_lim =  node->y_lim - pow_2[node->left_down->depth];
  node->left_down-> width = pow_2[node->depth];
  node->left_down-> cm_x = 0;
  node->left_down-> cm_y = 0;
  node->left_down-> tot_mass = -1.0;
  node->left_down->left_down = NULL;
  node->left_down->left_up = NULL;
  node->left_down->right_down = NULL;
  node->left_down->right_up = NULL;

}
//modify what node we are using (ie pointing to)
void insert(node_t **node, double x_pos, double y_pos, double mass, double* pow_2, int id){
  // inserts given bodyin tree, creating new nodes if necessary such that
  //every external node (leaf) contains only one body.
  //each node represents a quadrant and contains coordinates of
  //its center of mass and total mass within the quadrant.
  // N

  if((*node) == NULL){ printf("ERROR: Given node is NULL"); return;}

  //if node X does not contain a body, we put the new body in it.
  if((*node)->tot_mass == -1.0){
    //the node contains no bodies so we can insert the new one in it.
    (*node)-> cm_x = x_pos;
    (*node)-> cm_y = y_pos;
    (*node)-> tot_mass = mass;
    (*node)->body_id= id;

  }else{
    /*The node contains a body.
    If node is an internal node i.e has children we:
     update the center of mass and total mass of the node
     and go deeper to appropiate child (RU, RD, LU or LD) */
     if((*node)->right_up != NULL){
      //Node is internal (not a leaf)
      //UPDATE cm and MASS
      double cm_mass = (*node)->tot_mass;
      //update center of mass new_cm = (old_mass*old_cm + pos_new_particle*mass_new_particle)/(new total mass)
      (*node)->tot_mass+=mass;
      //we can save operations in here by declaring new varaibles
      (*node)->cm_x = cm_mass/(*node)->tot_mass * (*node)->cm_x + x_pos*(mass/(*node)->tot_mass); // Note: This maybe we can do more efficently
      (*node)->cm_y = cm_mass/(*node)->tot_mass* (*node)->cm_y + y_pos*(mass/(*node)->tot_mass);

      //FIND APPROPIATE CHILD and call insert with child as new node
      int expr = 2*(x_pos>(*node)->x_lim) + (y_pos >(*node)->y_lim);
      switch (expr) {
        case 0:
        //False, False we go left down
          node= &(*node)->left_down;
          break;
        case 1:
        //False, True we go left_up
          node= &(*node)->left_up;
          break;
        case 2:
        //True False, we go right down
          node = &(*node)->right_down;
          break;
        case 3:
        //True True we go right up
          node = &(*node)->right_up;
          break;
      }
      //insert with correct node
      insert(node, x_pos, y_pos, mass, pow_2, id);
    }else{
    //The node has mass but no children hence its a leaf
    //If it is an external node (leafs) we have
    //to create new children (furhter subdivide the space)
      create_children((*node), pow_2);
      //printf("We created children\n");
      //print_qtree((*node));
      //insert particle that was occupying the leaf, it will go to appropiate
      //quadrant
      double mass_in_node = (*node)->tot_mass;
      //set the mass to 0 since we remove a particle from the node
      (*node)->tot_mass = 0;
      //save position in the tree
      node_t ** parent_node = node;
      //insert particle that was already in node
      insert(parent_node, (*node)->cm_x, (*node)->cm_y, mass_in_node, pow_2,
      (*node)->body_id);
      //printf("We pushed particle that was already on the node\n");
      (*node)->body_id = -1;
      //print_qtree((*node));
      //try inserting particle i again
      insert(node, x_pos, y_pos, mass, pow_2, id);
      //printf("We inserted particle i \n");
      //print_qtree((*node));
      }
    }
    return;
  }
void delete_tree(node_t** node){
  if((*node)== NULL){return;}
  //If you can go to a child, ie child not null
  //go to it and call the function again
  if((*node)->right_up){
    delete_tree(&(*node)->right_up);
  }
  if((*node)->right_down){
    delete_tree(&(*node)->right_down);
  }
  if((*node)->left_up){
    delete_tree(&(*node)->left_up);
  }
  if((*node)->left_down){
    delete_tree(&(*node)->left_down);
  }
  //everything is NULL so we can free the node
  free((*node)->right_up);
  free((*node)->right_down);
  free((*node)->left_up);
  free((*node)->left_down);
  free((*node));
  (*node)=NULL;
  return;
}

void mark_as_not_use_tree(node_t** node){
    if((*node)== NULL || !(*node)->is_used){return;}
    //If you can go to a child, ie child not null
    //go to it and call the function again
    if((*node)->right_up && (*node)->right_up->is_used){
      //tree has a kid which is marked as used, so we go to it
      mark_as_not_use_tree(&(*node)->right_up);
    }
    if((*node)->right_down && (*node)->right_down->is_used){
      mark_as_not_use_tree(&(*node)->right_down);
    }
    if((*node)->left_up && (*node)->left_up->is_used){
      mark_as_not_use_tree(&(*node)->left_up);
    }
    if((*node)->left_down && (*node)->left_down->is_used){
      mark_as_not_use_tree(&(*node)->left_down);
    }
    //everything is NULL or marked as not used so we can mark this node
    //as unused
    (*node)->is_used = 0;
    return;
  }

typedef struct thread_arg{
  int thid;
  int i1;
  int i2;
  double theta_max;
  double G;
  double delta_t;
  double* pos_x;
  double* pos_y;
  double* vx;
  double* vy;
  node_t* node;
}arg_t;

void *  worker_get_acc(void* arg){
  /*
  This takes as arguments the root node of the tree,
  vector for positions and vector for velocities.
  Performs one update step for the particles i1 to i2-1 then deletes the tree.
  */
  arg_t* myarg = (arg_t*)arg;
  int i1 = myarg->i1;
  int i2 = myarg->i2;
  node_t* root = myarg->node;
  node_t** node;
  double G= myarg->G;
  double theta_max =myarg->theta_max;
  double total_acc_x=0;
  double total_acc_y=0;
  //printf("We are in worker for thread %d\n", myarg->thid);
  for(int i=i1; i <i2; i ++){
    node = &root;
    //calculate forces
    int c=0;
    while(root->is_used == 0){
      c++;
      //printf("We entered on while %d times \n",c );
      //printf("CURRENTLY ON NODE: %d, %lf, %lf \n", (*node)->depth, (*node)->cm_x,
      //(*node)->cm_y);
      //print_qtree(root);
      // Look if we are on a external node.
      double x_direction = (myarg->pos_x)[i] - (*node)->cm_x;
      //printf("x_direction : %lf \n", x_direction);
      double y_direction = (myarg->pos_y)[i] - (*node)->cm_y;
      //printf("y_direction: %lf \n", y_direction);
      double dist_to_node = sqrt(x_direction*x_direction + y_direction*y_direction);
      double tol= 1.0e-10; //
      if((*node)->is_used || dist_to_node < tol || (*node)->tot_mass < 0){
        //its the same particle or its already used
        (*node)->is_used = 1;
        node = &root;
      }else if((*node)->right_up &&
      (*node)->width/dist_to_node > theta_max){
        /*its internal node, and we need to go deeper*/
        int exp= (*node)->right_up->is_used + (*node)->right_down->is_used
        + (*node)->left_up->is_used + (*node)->left_down->is_used;
        //printf("exp: %d \n", exp);
        switch (exp) {
          case 0:
            //no child is used so we start with right up
            node = &(*node)->right_up;
            break;
          case 1:
            node = &(*node)->right_down;
            break;
          case 2:
            node = &(*node)->left_up;
            break;
          case 3:
            node = &(*node)->left_down;
            break;
          case 4:
          {
            //every child is used so we mark current node as used and go to root
            (*node)->is_used = 1;
            node = &root;
            break;
          } //we need brackets to define scope since we are declaring a variable
          default:
            printf("Something went quite wrong");
      }
    }else{
      double denominator = (dist_to_node+epsilon_0)*(dist_to_node+epsilon_0)*(dist_to_node+epsilon_0);
      //add to the global varaibles
      total_acc_x += G* (*node)->tot_mass*x_direction/denominator;
      total_acc_y += G* (*node)->tot_mass*y_direction/denominator;
      (*node)->is_used = 1;
      node = &root;
    }
  }
  //printf("Updating postion for particle %d, in thread %d \n", i, myarg->thid);
  //printf("Acceleration on %d: (%lf , %lf ) \n", i+1, total_acc_x, total_acc_y);
  myarg->vx[i] += myarg->delta_t  * total_acc_x;
  myarg->vy[i] += myarg->delta_t  * total_acc_y;
  myarg->pos_x[i] += myarg->delta_t * myarg->vx[i];
  myarg->pos_y[i] += myarg->delta_t * myarg->vy[i];
  //printf("pos: (%lf, %lf) \n", myarg->pos_x[i], myarg->pos_y[i]);
  //set acceleration to 0 and mark the tree as not used
  //once again for the next particle
  total_acc_x=0;
  total_acc_y=0;
  mark_as_not_use_tree(&root);
}
//we are done with the tree so we delete it.
delete_tree(&root);
//pthread_exit((void*) (intptr_t) myarg->thid);
}

int main(int argc, char *args[]){
  if (argc!=8){
    printf("Invalid number of arguments \n");
    printf("galsim expects: N filename nsteps delta_t theta_max graphics NUM_THREADS\n");
    return -1;
  }
  clock_t begin = clock();
  char *file_name = args[2];
  const int N = atoi(args[1]);
  const int n_steps = atoi(args[3]);
  /*not sure if this is the correct way of converting from
  character to double, maybe a single cast would suffice */
  const double delta_t = atof(args[4]);
  const double theta_max = atof(args[5]);
  const int NUM_THREADS = atoi(args[7]);
  const double G = -100/(double)N;
  //Read the file with initial conditions
  FILE *file;
  file = fopen(file_name , "rb");
  /*maybe in this case we could allocate memory for this
  matrix statically*/
  if(file_name == NULL){
    printf("ERROR reading the file");
  }

  /*maybe in this case we could allocate memory for this
  matrix statically*/
  double *pos_x = (double *)malloc(N*sizeof(double));
  double *pos_y = (double *)malloc(N*sizeof(double));
  double *vel_x = (double *)malloc(N*sizeof(double));
  double *vel_y = (double *)malloc(N*sizeof(double));
  double *ma = (double *)malloc(N*sizeof(double));
  double *bri = (double *)malloc(N*sizeof(double));

  for (int i = 0 ; i<(N) ; i++){
    double x , y , vx , vy , mass , bright;
    fread(&x , sizeof(double) , 1 ,file);
    fread(&y , sizeof(double) , 1 ,file);
    fread(&mass , sizeof(double) , 1 ,file);
    fread(&vx , sizeof(double) , 1 ,file);
    fread(&vy , sizeof(double) , 1 ,file);
    fread(&bright, sizeof(double) , 1 ,file);
    pos_x[i] = x;
    pos_y[i] = y;
    ma[i] = mass;
    vel_x[i] = vx;
    vel_y[i] = vy;
    bri[i] = bright;
  }
  fclose(file);
  printf("We read the file \n");
  //find powers of two so we only have to do it once
  const int K=200;
  double pow_2[K];
  pow_2[0]=1.0;
  for(int i=1; i< K; i++)
    pow_2[i]= pow_2[i-1]/2;

  //initialize root
  // insert particles in quad tree one by one
  //we pass a pointer to the node since we want to
  //modify what node we are using (ie pointing to)

  //printf("Particles: \n");
  //for(int i=0; i<N; i++){
  //  printf("particle %d: m = %lf, (%lf,%lf) \n", i, arr[i][MASS],
  //  arr[i][POS_X], arr[i][POS_Y]);
  //}
  //printf("Printing tree .. \n");
  //print_qtree(root);

  for (int k = 0 ; k<n_steps ; k++){
    //printf("We are on step %d\n", k);
    //create as many trees as threads
    //printf("Allocating mem for root \n");
    node_t** trees_root = (node_t**)malloc(sizeof(node_t*)*NUM_THREADS);
    for(int i=0; i< NUM_THREADS; i++){
      node_t* root = (node_t*)malloc(sizeof(node_t));
      root->depth = 1;
      root->body_id = i;
      root->x_lim = 0.5;
      root->y_lim = 0.5;
      root->is_used = 0;
      root->width = 0.5;
      root->tot_mass = -1.0;
      root->left_down = NULL;
      root->left_up = NULL;
      root->right_up=NULL;
      root->right_down=NULL;
      trees_root[i] = root;
      //printf("Created root %d ... \n", i);

    /*  int depth;
      int body_id;
      int is_used;
      double x_lim; //x division (middle point in x)
      double y_lim; // y division (middle point in y)
      double width; //width of the box
      double cm_x; // = 0, this points to 0 and
      double cm_y; //center of mas s of the quadrant
      double tot_mass; // = 0, total mass of the quadrant
      struct tree_node *left_down; //Q3 child
      struct tree_node *left_up; //Q2 child
      struct tree_node *right_down; //Q4 child
      struct tree_node *right_up; //Q1 child*/
    }
    /// CHECKING IF ITS OK
    arg_t* arg_thread = (arg_t*)malloc(NUM_THREADS*sizeof(arg_t));
    for(int t = 0; t< NUM_THREADS; t++){
      printf("Root %d, depth: %d, %lf. Pls be different: %d\n", t, trees_root[t]->depth,
      trees_root[t]->x_lim, trees_root[t]->body_id );
    }
    int id_tar;
    //Create the trees for step K
    for(int j =0; j < NUM_THREADS; j++){
      for(int i=0; i<N; i++){
        id_tar = i;
        insert(&(trees_root[j]),pos_x[i], pos_y[i], ma[i], pow_2, id_tar);
      }
    }
    //printf("Inserted nodes ...\n");
    //Start threads
    int rc;
    int t;

    pthread_t *thread = (pthread_t*)malloc(sizeof(pthread_t)*NUM_THREADS);

    pthread_attr_t attr;
    /* Initialize thread attr and set to JOINABLE*/
    pthread_attr_init(&attr); //initializes thread attributes with default values
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE); //set it to detached

    //printf("ceating threads \n");
    for(t=0; t<NUM_THREADS; t++) {
       printf("Main: creating thread %d\n", t);
       //prepare thread arguments
       arg_t* thread_arg = &arg_thread[t];
       thread_arg->thid = t;
       thread_arg->i1 = N/NUM_THREADS*t;
       thread_arg->i2 = N/NUM_THREADS*(t+1);
       thread_arg->theta_max = theta_max;
       thread_arg->G = G;
       thread_arg->delta_t = delta_t;
       thread_arg->pos_x = &pos_x[0];
       thread_arg->pos_y = &pos_y[0];
       thread_arg->vx = &vel_x[0];
       thread_arg->vy= &vel_y[0];
       thread_arg->node = trees_root[t];
       rc = pthread_create(&thread[t], &attr, worker_get_acc, thread_arg);
       /*if (rc) {
        printf("ERROR; return code from pthread_create() is %d\n", rc);
        exit(-1);
      }*/

    }
    //free(arg_thread);


    // MAIN CAN TAKE CARE OF THE REST. (calling worker again)
    //destroy thread attr and join before deleting tree
    pthread_attr_destroy(&attr);
    void* status;
    for(int t=0; t<NUM_THREADS; t++) {
       rc = pthread_join(thread[t], &status);
       if (rc) {
          //printf("ERROR; return code from pthread_join() is %d\n", rc);
          exit(-1);
          }
      // printf("Main: completed join with thread %d having a status of %ld\n",t,(long)status);
      }
    //delete trees THIS IS probably better to put it on the worker function for each
    //thread
    free(thread);
    free(trees_root);
    free(arg_thread);
  }

  FILE * fout = fopen("result.gal", "w+");          //check succesful creation/opening of results file
    if (fout == NULL){
      printf("ERROR opening the output file");
      exit(1);                                    // added exit to stop prgram too
    }

// print_struct(n,galaxy);

for (int j=0; j<N; j++){
  double x , y , vx , vy;
  x = pos_x[j];
  y = pos_y[j];
  vx = vel_x[j];
  vy = vel_y[j];
  fwrite(&x, sizeof(double), 1, fout);
  fwrite(&y,  sizeof(double), 1, fout);
  fwrite(&ma[j] ,  sizeof(double), 1, fout);
  fwrite(&vx , sizeof(double), 1, fout);
  fwrite(&vy,  sizeof(double), 1, fout);
  fwrite(&bri[j],  sizeof(double), 1, fout);
}

// print_file(N, fout);
// print_file(N,fout);
fclose(fout);

  free(pos_x);
  free(pos_y);
  free(vel_x);
  free(vel_y);
  free(ma);
  free(bri);


  return 0;
}
