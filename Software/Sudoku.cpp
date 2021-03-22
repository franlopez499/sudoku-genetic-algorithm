
/* ------------------------- PROBLEMA DE LAS N REINAS ----------------------- */
#include <ga/GASimpleGA.h> //  Algoritmo Genetico simple
#include <ga/GA1DArrayGenome.h> // Genoma --> array de enteros (dim. 1) alelos
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;


#define MAX_GENERACIONES 12000
struct plantilla{
       int tam;   
       int *fijo;   
};

float Objective(GAGenome &); // Funcion objetivo 

GABoolean Termina(GAGeneticAlgorithm &); // Funcion de terminacion 

void InicioSudoku(GAGenome& g);
void leerSudoku(struct plantilla *S,char *nombreF);
bool checkColumna(int col[], int * check, int tam);
int MutacionSudoku(GAGenome& g,float pmut);
int CruceSudoku(const GAGenome& p1,const GAGenome & p2,GAGenome* c1,GAGenome* c2);


int main(int argc, char **argv)
{

    char * selector = argv[1];
    int popsize = atoi(argv[2]);
    float pcross = atof(argv[3]);
    float pmut = atof(argv[4]);
    char * filepath = argv[5];

    struct plantilla sudoku;
    leerSudoku(&sudoku, filepath);

    // Conjunto enumerado de alelos --> valores posibles de cada gen del genoma

    GAAlleleSet<int> alelos;
    for(int i=1;i <= sudoku.tam; i++) alelos.add(i);

    // Creamos el genoma y definimos operadores de inicio, cruce y mutaci�n

    GA1DArrayAlleleGenome<int> genome(sudoku.tam*sudoku.tam,alelos,Objective,&sudoku);
    genome.crossover(CruceSudoku);
    genome.mutator(MutacionSudoku);
    genome.initializer(InicioSudoku);

// Creamos el algoritmo genetico

    GASimpleGA ga(genome);

    // Inicializamos - minimizar funcion objetivo, tama�o poblacion, n� generaciones,
// pr. cruce y pr. mutacion, selecci�n y le indicamos que evolucione.
    ga.minimaxi(-1);
    ga.populationSize(popsize);
    ga.nGenerations(MAX_GENERACIONES);
    ga.pCrossover(pcross);
    ga.pMutation(pmut);
    bool s;
    GARouletteWheelSelector roulette_selector;
    GATournamentSelector tournament_selector;

    if(!strcmp("-r",selector)) {ga.selector(roulette_selector); s = true;}
    else if(!strcmp("-t",selector)){ ga.selector(tournament_selector); s = false;}
    else {
        fprintf(stderr,"Error: selector '%s' no reconocido\n",selector);
        exit(EXIT_FAILURE);
    }

    ga.terminator(Termina);
    ga.evolve(1);
    cout <<"\t"<< (s ? "Selector: GARouletteWheelSelector" : "Selector: GATournamentSelector")<< endl
    << "\tTam de poblacion: " << popsize << endl << "\tProbabilidad de cruce: " << pcross << endl << "\tProbabilidad de mutacion: "
    << pmut << endl << "\tArchivo: " << filepath << endl<< "\tFitness: "<<ga.statistics().minEver() << endl
    <<  "\tNumero de generaciones: " << ga.statistics().generation() << endl <<"\t"<<(ga.statistics().minEver()==0?"Estado: RESUELTO\t" : "Estado: NO RESUELTO\t") <<endl<< endl;
    //cout << (s ? "GARouletteWheelSelector" : "GATournamentSelector") << ";" << popsize << ";" << pcross << ";" << pmut  << ";"<<ga.statistics().minEver() <<  ";" << ga.statistics().generation() << endl;

    cout << "Estado final del sudoku:\n\n";
    GA1DArrayAlleleGenome<int> & best = (GA1DArrayAlleleGenome<int> &)ga.statistics().bestIndividual();
    for(int i = 0; i<sudoku.tam ; i++){
        cout << "\t[ ";
        for(int j=0; j< sudoku.tam; j++){
            cout << best.gene(i*sudoku.tam+j) ;
            cout << " ";
        }
        cout << "]"<<endl;

    }
    cout << endl<< "---------------------------------"<<endl;


}



// Funcion de terminacion

GABoolean Termina(GAGeneticAlgorithm & ga){
    if ( (ga.statistics().minEver()==0) ||
        (ga.statistics().generation()==ga.nGenerations()) ) return gaTrue;
    else return gaFalse;
}

void leerSudoku(struct plantilla *S,char *nombreF){
   ifstream f(nombreF);

   f>>S->tam;
   S->fijo = new int[S->tam*S->tam];

   for(int i=0;i<S->tam*S->tam;i++){
        f>>S->fijo[i];
    }

   f.close();
}

void InicioSudoku(GAGenome& g){


     //Casting al individuo
     GA1DArrayAlleleGenome<int> & genome = ( GA1DArrayAlleleGenome<int> &)g;

     struct plantilla * plantilla1;
     plantilla1 = (struct plantilla *) genome.userData();

     int aux[plantilla1->tam];

     for(int f=0;f<plantilla1->tam;f++){

          for(int j=0;j<plantilla1->tam;j++) aux[j]=0;  // Inicializacion de fila a 0

          // Con este bucle inicializamos los numeros de una fila aleatoriamente sin repeticiones
          for(int j=1;j<=plantilla1->tam;j++){
            int v=GARandomInt(0,plantilla1->tam-1);
            while (aux[v]!=0) v=(v+1)%plantilla1->tam; // Buscamos una posicion libre
            aux[v]=j;   //Almacenamos el valor en la posicion correspondiente
          }

          int i=0;

          while(i<plantilla1->tam){

              while((plantilla1->fijo[(f*plantilla1->tam)+i]==0) && (i<plantilla1->tam)) i++;

              if (i<plantilla1->tam){

                     bool encontrado=false;
                     for(int j=0;(j<plantilla1->tam) && (!encontrado);j++)
                             if (aux[j]==plantilla1->fijo[(f*plantilla1->tam)+i]) {
                                encontrado=true;
                                aux[j]=aux[i];  // Intercambiamos valores
                             }

                     aux[i]=plantilla1->fijo[(f*plantilla1->tam)+i];
              }
              i++;

          }

          for(int c=0;c<plantilla1->tam;c++)
            genome.gene((f*plantilla1->tam)+c,aux[c]);
     }
}

int CruceSudoku(const GAGenome& p1,const GAGenome & p2,GAGenome* c1,GAGenome* c2){


    // Casting
    GA1DArrayAlleleGenome<int> & m = ( GA1DArrayAlleleGenome<int> &)p1;
    GA1DArrayAlleleGenome<int> & p = ( GA1DArrayAlleleGenome<int> &)p2;

    struct plantilla * plantilla1 = (struct plantilla *) m.userData();
    int n=0;

    int punto1=GARandomInt(0,m.length());   //Genera un punto de corte aleatorio dentro del tama�o del sudoku
    while ((punto1%plantilla1->tam)!=0) punto1++; // Tiene que ser modulo 0 para asi poder copiar desde el principio de la fila, respetando la no repeticion de elementos en filas
    int punto2=m.length()-punto1;   //Obtenemos en punto2 el tama�o que queda por recorrer, el cual usaremos en la funcion copy para copiar los genes restantes del padre en hx

    if (c1){


            // Casting a GA1DArrayGenome ya que copy() es miembro de esta clase
            GA1DArrayGenome<int> &h1 = (GA1DArrayGenome<int> &)*c1;

             h1.copy(m,0,0,punto1);   // Se copian al hijo 1 los primeros punto1 elementos de la madre
             h1.copy(p,punto1,punto1,punto2); //Copiamos al hijo 1 desde el elemento punto1 hasta el final del padre.
             n++;
    }

    if (c2){

            // Casting a GA1DArrayGenome ya que copy() es miembro de esta clase
            GA1DArrayGenome<int> &h2=(GA1DArrayGenome<int> &)*c2;

             h2.copy(p,0,0,punto1);
             h2.copy(m,punto1,punto1,punto2);
             n++;
    }

    return n;

}

bool checkColumna(int col[], int * check, int tam){
     bool repe=false;

    for(int i=0;i<tam;i++) check[i]=0;

    for(int i=0;i<tam;i++)
             check[col[i]-1]++;

    for(int i=0;i<tam;i++) if (check[i]>1) repe=true;

     return repe;
}


int MutacionSudoku(GAGenome& g,float pmut){

    // Casting
    GA1DArrayAlleleGenome<int> & genome = ( GA1DArrayAlleleGenome<int> &)g;

    struct plantilla * plantilla1;
    plantilla1 = (struct plantilla *) genome.userData();
    int nmut=0;
    int aux;
    int fil = 0;
    bool fila;

    int caux[plantilla1->tam];
    int *checkC=new int[plantilla1->tam];
    if (pmut<=0.0) return 0;

    for(int f=0; f<plantilla1->tam; f++)
       for(int c=0; c<plantilla1->tam; c++)
          if (plantilla1->fijo[(f*plantilla1->tam)+c]==0){
           if (GAFlipCoin(pmut) ){
                if (GAFlipCoin(0.5)) fila = true;
                else fila = false;

                if (!fila){ //Mutacion por columnas

                      for(int j=0;j<plantilla1->tam;j++) caux[j]=genome.gene((j*plantilla1->tam)+c);
                      if (checkColumna(caux,checkC,plantilla1->tam)){
                         int v1 = GARandomInt(0,plantilla1->tam-1);
                         while (checkC[v1]<=1) v1=(v1+1)%plantilla1->tam; //Buscamos un elemento repetido
                         v1++;
                         int v2 = GARandomInt(0,plantilla1->tam-1);
                         while (checkC[v2]!=0) v2=(v2+1)%plantilla1->tam; //Buscamos un elemento no repetido
                         v2++;

                         bool encontrado = false;
                         for(int j=0;j<plantilla1->tam && !encontrado;j++)
                                 if ((plantilla1->fijo[j*(plantilla1->tam)+c]==0)&&(genome.gene(j*(plantilla1->tam)+c)==v1)){ // Comprobamos que no sea un dato inicial y que esa posicion contenga el valor v1
                                    encontrado = true;
                                    genome.gene((j*plantilla1->tam)+c,v2);    //Coloco v2, elemento no repetido en esa posicion
                                    fil = j;
                                 }

                         int col=(c+1)%plantilla1->tam;
                         while(genome.gene((fil*plantilla1->tam)+col)!=v2) col=(col+1)%plantilla1->tam;   // Buscamos la posicion donde v2 se repite
                         if (plantilla1->fijo[(fil*plantilla1->tam)+col]==0) { //Si la pos de la plantilla en la que estamos es variable
                                nmut++;
                                genome.gene((fil*plantilla1->tam)+col,v1);  // Intercambiamos el valor repetido por v1 para asegurarnos de que la mutacion no genera repeticiones en la fila
                         }
                         else {                                          //Si es un dato inicial restauramos el valor para que no haya repeticiones en las filas
                              genome.gene((fil*plantilla1->tam)+c,v1);
                         }

                      }

                }
                else{   //Mutacion por filas
                   int v1 = (c + 1) %plantilla1->tam;
                   while ((plantilla1->fijo[(f*plantilla1->tam)+v1]!=0)) v1=(v1+1)%plantilla1->tam;
                   aux = genome.gene((f*plantilla1->tam)+c);
                   genome.gene((f*plantilla1->tam)+c,genome.gene((f*plantilla1->tam)+v1));
                   genome.gene((f*plantilla1->tam)+v1,aux);
                   nmut++;
                }
           }
          }

    return nmut;
}

float Objective(GAGenome& g) {
    GA1DArrayAlleleGenome<int> & genome = (GA1DArrayAlleleGenome<int> &)g;
    struct plantilla * sudoku;
    sudoku = (struct plantilla *) genome.userData();
    int tam = sudoku->tam;
    float fitness = 0.0f;


    //Vector que almacena el numero de repeticiones de cada numero
    vector<int> conjunto(tam,0);


    for(int j = 0; j<tam ; j++)
    {
        for(int i=0; i< tam; i++)
        {
            int elemento = genome.gene(i*tam+j)-1;
            if(conjunto[elemento]++)
                fitness+=1.0f;
        }
        fill(conjunto.begin(),conjunto.end(),0);
    }


    //Contar repeticiones en cuadrantes
    int tamCuandrante = sqrt(tam);
    for(int x =0; x<tamCuandrante; x++)
    {
        for(int y =0; y<tamCuandrante; y++)
        {
            for(int i = 0; i< tamCuandrante ; i++)
            {
                for(int j=0; j< tamCuandrante; j++)
                {
                    int fila = i+x*tamCuandrante;
                    int columna = j+y*tamCuandrante; // SOLO AVANZA I E J
                    int elemento = genome.gene(fila*tam+columna)-1;
                    if(conjunto[elemento]++)
                        fitness+=1.0f;

                }
            }
            fill(conjunto.begin(),conjunto.end(),0);
        }
    }


    return fitness;
}
