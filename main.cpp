/*
*  main.cpp
*
*  M�dulo pricipal do estudo TSP
*  Autor: Geraldo Jos� Ferreira Chagas Junior - gjr.doc@gmail.com
*
*  PPGI - NCE - UFRJ
*  Data Cria��o: 05/05/2017
*  Datas de Modifica��es:
*
* Programa de Algor�tmo Genetico para resolver tsp de forma recursiva
* Necess�rio instala��o da biblioteca libxml
* Instala��o da biblioteca no debian
*           $> sudo apt-get install libxml
*
*  Os parametros de entrada est�o explicados no arquivo tsphelp.txt
*  A compila��o pode ser realizada pelo comando make
*
*  Se n�o for definido parametros de entrada, ser�o utilizados os
*  parametros padr�es, conforme abaixo.
*
*  n�mero de execu��es:
*  muta��o:
*  cruzamento:
*  popula��o:
*  % de muta��o:
*  % de cruzamento:
*  recursividade:
*  % de melhora recursiva:
*  arquivo de saida:
*
*/
#include <iostream>
#include <string>
#include "config.hpp"
#include "ag.hpp"
#include "tsp.hpp"

#define strXML "-XML"
#define strTSP "-TSP"

using namespace std;

#ifdef LIBXML_TREE_ENABLED

int main(int argc, char *argv[])
{
   string nomeArqSaida;
   string cabecalho;
   string typeArq;

   locale loc;

   TAlgGenetico *ag;
   TMapaGenes *mapa = new TMapaGenes();
   TConfig *config  = new TConfig();
   TArqLog *arqSaida;

   time_t tempo;
   struct tm *tlocal;
   char data[128];

   //par�metros obrigat�ros como entrada
   if (argc < 5)
   {
      cout << "Par�metros obrigat�ros:" << endl;
      cout << "\t 1 - -XML ou -TSP. (arquivo de inst�ncia XML ou TSP)" << endl;
      cout << "\t 2 - Arquivo de inst�ncia TSP, no formato XML OU TSP" << endl;
      cout << "\t 3 - Arquivo de configura��o, no formato XML" << endl;
      cout << "\t 4 - Nome do arquivo de sa�da, resultados" << endl;
      return 1;
   }

   LIBXML_TEST_VERSION

   typeArq = argv[1];

   cout << "Caregando arquivo de configura��o: " << argv[3] << endl;
   config->carregaDoArquivo(argv[3]);
   cout << "Arquivo " << argv[3] << " carregado." << endl;

   cout << "Caregando inst�ncia " << argv[2] << endl;

   if (0 == typeArq.compare(strXML))
   {
      mapa->carregaDoArquivo (argv[2]);
   }
   else if ((0 == typeArq.compare(strTSP)))
   {
       if (!mapa->carregaDoArquivoTSP (argv[2]))
          return 1;
   }
   else
   {
      cout << "Formato de arquivo desconhecido " << typeArq << endl;
      return 1;
   }

   cout << "Inst�ncia " << argv[2] << " carregada." << endl;

   TUtils::initRnd ();

   for (int countExec=0; countExec<config->numExec; countExec++)
   {
      /* para complemento do nome do arquivo de saida */
      tempo = time(0);
      tlocal = localtime(&tempo);
      strftime(data, 128, "%d_%m_%y_%H_%M_%S", tlocal);
      nomeArqSaida = argv[4];
      nomeArqSaida += "_";
      nomeArqSaida += to_string(countExec);
      nomeArqSaida += "_";
      nomeArqSaida += data;
      nomeArqSaida += ".txt";

      if (config->printParcial)
      {
         cout << "Execu��o " << countExec+1 << " / " << config->numExec << " iniciada. "<< endl;
         cout << "Arquivo de saida de resultados: " << nomeArqSaida << endl;
      }

      cabecalho   = "Inst�ncia;";
      cabecalho  +=  argv[2];
      cabecalho  +=  "\n";

      cabecalho  += "Execu��o;";
      cabecalho  += to_string(countExec+1);
      cabecalho  += " / ";
      cabecalho  += to_string(config->numExec);
      cabecalho  += "\n";

      cabecalho  += "Tamanho da Popula�ap;";
      cabecalho  += to_string(config->tamPopulacao);
      cabecalho  += "\n";

      cabecalho  += "M�ximo de gera��es;";
      cabecalho  += to_string(config->maxGeracao);
      cabecalho  += "\n";

      cabecalho  += "Muta��o;";
      cabecalho  += to_string(config->mutacao);
      cabecalho  += "\n";

      cabecalho  += "Cruzamento;";
      cabecalho  += to_string(config->cruzamento);
      cabecalho  += "\n";

      cabecalho  += "% Manpula��o;";
      cabecalho  += to_string(config->percentManipulacao);
      cabecalho  += "\n";

      cabecalho  += "% Muta��o;";
      cabecalho  += to_string(config->percentMutacao);
      cabecalho  += "\n";

      cabecalho  += "Print Parcial;";
      cabecalho  += to_string(config->printParcial);
      cabecalho  += "\n";

      cabecalho  += "Ativa Recursivo;";
      cabecalho  += to_string(config->percentMutacaoRecursiva);
      cabecalho  += "\n";

      cabecalho  += "Percentual de Redu��o;";
      cabecalho  += to_string(config->percentReducao);
      cabecalho  += "\n";

      cabecalho  += "Profundidade M�xima;";
      cabecalho  += to_string(config->profundidadeMaxima);
      cabecalho  += "\n";

      cabecalho  += "Percentual de Eltismo;";
      cabecalho  += to_string(config->percentElitismo);
      cabecalho  += "\n";

      cabecalho  += "Sele��o para Cruzamento;";
      cabecalho  += to_string(config->selecao);
      cabecalho  += "\n";

      cabecalho  += "Forma de Sele��o par Muta��o;";
      cabecalho  += to_string(config->selIndMutacao);
      cabecalho  += "\n";

      arqSaida = new TArqLog(cabecalho, nomeArqSaida);
      ag = new TAlgGenetico(mapa, arqSaida);
      ag->setMutacao(config->mutacao);
      ag->setCruzamento(config->cruzamento);
      ag->setTamPopulacao(config->tamPopulacao);
      ag->setPrintParcial(config->printParcial);
      ag->setMaxGeracao(config->maxGeracao);
      ag->setPercentElitismo(config->percentElitismo);
      ag->setPercentMutacao(config->percentMutacao);
	   ag->setProfundidadeMaxima(config->profundidadeMaxima);
      ag->setSelecao(config->selecao);
      ag->setSelIndMutacao(config->selIndMutacao);
      ag->setPercentMutacaoRecursiva(config->percentMutacaoRecursiva);
      ag->setPercentReducao(config->percentReducao);
      ag->setTimeOut(config->timeout);
      ag->exec();

      arqSaida->addLinha("");
      arqSaida->addLinha("");

      delete arqSaida;
      delete ag;
   }

   delete config;
   delete mapa;
   return 0;
}
#else
int main(void)
{
    cout << "Suporte a �rvore xml n�o compilado" << endl;
    exit(1);
}
#endif
