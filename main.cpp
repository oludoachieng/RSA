#include <stdio.h>
#include <gmp.h>
#include <gmpxx.h>
#include <assert.h>
#include<cstdlib>
#include <fstream>
#include <iostream>
#include<ctime>
#include <string>
#include <cstring>
#include <sstream>
using namespace std;

int main(){
//Prompt user for either of the the choices
  int selection,counter=0, trial=0;
  printf ("Press 1 to generate parameters. ");
  printf ("\n");
  printf ("Press 2 to encrypt. ");
  printf ("\n");
  printf ("Press 3 to decrypt. ");
  printf ("\n");
  scanf("%d" , &selection); //take input value

//generate parameters
if(selection==1){

new_rand:
//generate random 512 bit number
mpz_class gen_no,p,q;
gmp_randclass ran(gmp_randinit_default);
ran.seed(time(NULL)+rand());
gen_no =ran.get_z_bits(512);

//random number completed

//Ensure the number generated is odd
if(gen_no % 2==0)
goto new_rand;
//print out generated number


check_prime:
//Miller Rabin primality test (same as homework 1)to check if our number is prime
 // Set the value of d as n-1
  mpz_t d;
  mpz_init(d);
  mpz_set_ui(d,0);
  mpz_sub_ui(d,gen_no.get_mpz_t(),1);
 //mpz_out_str(stdout,10,d);

  //use the mpz class in order to use C++ primitives for loops and arithmetic
  mpz_class v (d);
  mpz_class s (0);

 //computing n-1=(2^s)*v where v is an odd number
  do
   {
     v=v/2;
     s=s+1;
   }
   while(v% 2==0);

   gen_rand:

//generate a random number 'a' mod d

srand(time(0));

more_trials: //Shall use goto loop to repeat test severally at this point

mpz_class a=rand();
mpz_class expo (1);
mpz_class expo2 (2);

//Reduce a mod n
mpz_powm (a.get_mpz_t(), a.get_mpz_t(), expo.get_mpz_t(), d);

//check if the random number generator gives zero
    if(a==0)
    goto gen_rand;

//Assign b=a^v mod n
mpz_class b;
mpz_powm (b.get_mpz_t(), a.get_mpz_t(), v.get_mpz_t(), gen_no.get_mpz_t());


//check if b%n=1
mpz_class rem_1;
mpz_powm (rem_1.get_mpz_t(), b.get_mpz_t(), expo.get_mpz_t(), gen_no.get_mpz_t());
if(rem_1==1)//probably prime
{
   goto trial_loop; //we have to repeat with different a's to reduce error
}

//for loop from 0 to s-1 where s is the exponent of 2
for(int i=0;i<s;i++)
{
mpz_class rem_2(b+1);
mpz_powm (rem_2.get_mpz_t(), rem_2.get_mpz_t(), expo.get_mpz_t(), gen_no.get_mpz_t());

    if(rem_2==0)//probably prime
      {
       goto trial_loop;//we have to repeat with different a's to reduce error
      }
    else//Assign b=b^2 mod n
        {
         mpz_powm (b.get_mpz_t(), b.get_mpz_t(), expo2.get_mpz_t(), gen_no.get_mpz_t());
        }
}
//If test fails increase our number by 2 and check if new number is prime
  mpz_add_ui(gen_no.get_mpz_t(),gen_no.get_mpz_t(),2);
  goto check_prime;


//For probably prime repeat test 5 times to reduce error to approximately(9.7*10^-4)
trial_loop:
if(trial<6)
{
    trial++;
    goto more_trials;
}

//Write the generated p value to parameters.txt
if(counter==0)
{
// open a file in write mode.
   ofstream outfile;
   outfile.open("parameters.txt");

   // write inputted data into the file.
   outfile << gen_no.get_mpz_t()<<endl ;
    p=gen_no;

   // close the opened file.
   outfile.close();
}

//Write the generated q value to parameters.txt

  counter++;
  while(counter<2)
  {
    goto new_rand;
     }
  ofstream outfile;
  outfile.open("parameters.txt", fstream::app);
  outfile << gen_no.get_mpz_t()<<endl ;
  q=gen_no;
  outfile.close();


ifstream infile;
string line,line2;
mpz_t num,num2;

//reading from the text file
infile.open("parameters.txt");
if(infile.is_open()){

   if (infile.good()){
       getline(infile, line);
       getline(infile, line2);

    }
}
infile.close();

char cstr[line.size() + 1];
	strcpy(cstr, line.c_str());

 //Initialize the number num
  mpz_init(num);
  mpz_set_ui(num,0);

//Parse the input string as a base 10 number
  mpz_set_str(num,cstr, 10);
  mpz_class w (num);

 char cstr2[line2.size() + 1];
	strcpy(cstr2, line2.c_str());

 //Initialize the number num2
  mpz_init(num2);
  mpz_set_ui(num2,0);

//Parse the input string as a base 10 number
  mpz_set_str(num2,cstr2, 10);
  mpz_class y (num2);

//n=p*q
mpz_class num3;
mpz_mul(num3.get_mpz_t(),w.get_mpz_t(),y.get_mpz_t());

//subtract p-1 and q-1
mpz_sub_ui(w.get_mpz_t(),w.get_mpz_t(),1);
mpz_sub_ui(y.get_mpz_t(),y.get_mpz_t(),1);

//multiply p-1 and q-1 to get phi
mpz_class phi;
mpz_mul(phi.get_mpz_t(),w.get_mpz_t(),y.get_mpz_t());
printf("\n");

//append phi to parameters.txt
  outfile.open("parameters.txt", fstream::app);
  outfile << phi.get_mpz_t()<<endl ;
  outfile.close();
//append n to parameters.txt
outfile.open("parameters.txt", fstream::app);
  outfile << num3.get_mpz_t()<<endl ;
  outfile.close();
//generate e whose gcd with phi is 1
regenerate:
mpz_class gen_e,gcd_e;
gmp_randclass ran_e(gmp_randinit_default);
ran_e.seed(time(NULL)+rand());
gen_e =ran.get_z_bits(512);

mpz_gcd (gcd_e.get_mpz_t(), phi.get_mpz_t(), gen_e.get_mpz_t());//get gcd( (phi),e)
if(gcd_e!=1){goto regenerate;}

outfile.open("parameters.txt", fstream::app);
  outfile << gen_e.get_mpz_t()<<endl ;
  outfile.close();
  //find d which is e inverse mod phi
  mpz_class dee;
  mpz_invert(dee.get_mpz_t(), gen_e.get_mpz_t(), phi.get_mpz_t());
outfile.open("parameters.txt", fstream::app);
  outfile << dee.get_mpz_t()<<endl ;
  outfile.close();

cout<<"Parameters generated and stored in parameters.txt file in order: p, q, phi, n, e and d."<<endl;

return 0;
}


//encrypt plaintext
else if(selection==2){
     int length,length2;
    ifstream filestr;
//Error if no plaintext and parameters
filestr.open("plain.txt", ios::binary); // open your file
filestr.seekg(0, ios::end); // put the "cursor" at the end of the file
length = filestr.tellg(); // find the position of the cursor
filestr.close(); // close your file

filestr.open("parameters.txt", ios::binary); // open your file
filestr.seekg(0, ios::end); // put the "cursor" at the end of the file
length2 = filestr.tellg(); // find the position of the cursor
filestr.close(); // close your file


if ( length == 0 )
{
    cout<<"Enter plaintext to plain.txt.";
    return 0;
    }
if ( length2 == 0 )
{
    cout<<"Generate parameters first.";
    return 0;
    }


else {
string encr;
mpz_t numb;

ifstream plainfile;
//reading from the text file
plainfile.open("plain.txt");
if(plainfile.is_open()){

   if (plainfile.good()){
       getline(plainfile, encr);
    }
}
plainfile.close();
//convert to mpz class
char enstr1[encr.size() + 1];
strcpy(enstr1, encr.c_str());
  mpz_init(numb);
  mpz_set_ui(numb,0);
  mpz_set_str(numb,enstr1, 10);
  mpz_class v (numb);

//capture all values in the parameters.txt that will be used
ifstream infile;
string line1,line2,line3,line4,line5,line6;
mpz_t num1,num2,num3,num4,num5,num6;


//reading from the text file
infile.open("parameters.txt");
if(infile.is_open()){

   if (infile.good()){
       getline(infile, line1);
       getline(infile, line2);
       getline(infile, line3);
       getline(infile, line4);
       getline(infile, line5);
       getline(infile, line6);
    }
}
infile.close();
//capture p
char mestr1[line1.size() + 1];
strcpy(mestr1, line1.c_str());
  mpz_init(num1);
  mpz_set_ui(num1,0);
  mpz_set_str(num1,mestr1, 10);
  mpz_class a (num1);
//capture q
  char mestr2[line2.size() + 1];
strcpy(mestr2, line2.c_str());
  mpz_init(num2);
  mpz_set_ui(num2,0);
  mpz_set_str(num2,mestr2, 10);
  mpz_class b(num2);
//capture phi
  char mestr3[line3.size() + 1];
strcpy(mestr3, line3.c_str());
  mpz_init(num3);
  mpz_set_ui(num3,0);
  mpz_set_str(num3,mestr3, 10);
  mpz_class c (num3);
//capture n
  char mestr4[line4.size() + 1];
strcpy(mestr4, line4.c_str());
  mpz_init(num4);
  mpz_set_ui(num4,0);
  mpz_set_str(num4,mestr4, 10);
  mpz_class d (num4);
//capture e
  char mestr5[line5.size() + 1];
strcpy(mestr5, line5.c_str());
  mpz_init(num5);
  mpz_set_ui(num5,0);
  mpz_set_str(num5,mestr5, 10);
  mpz_class e (num5);
//capture d
  char mestr6[line6.size() + 1];
strcpy(mestr6, line6.c_str());
  mpz_init(num6);
  mpz_set_ui(num6,0);
  mpz_set_str(num6,mestr6, 10);
  mpz_class f (num6);

//calculate the ciphertext
mpz_t cipher;
mpz_init(cipher);
mpz_set_ui(cipher,0);
mpz_class cip (cipher);
//c=m power e mod n
mpz_powm (cip.get_mpz_t(), v.get_mpz_t(), e.get_mpz_t(), d.get_mpz_t());

  ofstream cipoutfile;
  cipoutfile.open("cipher.txt");
  cipoutfile << cip.get_mpz_t()<<endl ;
  cipoutfile.close();
  cout<<"Encryption complete and stored in cipher.txt"<<endl;
 return 0;
}
}

//Decrypt the ciphertext
else if(selection==3){
    int length,length2;
    ifstream filestr;
//Error if cipher.txt and parameters.txt are empty
filestr.open("cipher.txt", ios::binary); // open your file
filestr.seekg(0, ios::end); // put the "cursor" at the end of the file
length = filestr.tellg(); // find the position of the cursor
filestr.close(); // close your file

filestr.open("parameters.txt", ios::binary); // open your file
filestr.seekg(0, ios::end); // put the "cursor" at the end of the file
length2 = filestr.tellg(); // find the position of the cursor
filestr.close(); // close your file


if ( length == 0 )
{
    cout<<"Enter ciphertext to cipher.txt.";
    return 0;
    }
if ( length2 == 0 )
{
    cout<<"Generate parameters first.";
    return 0;
    }
//decryption begins
else{
ifstream decfile;
string cipha;
mpz_t decry;

//reading from the text file
decfile.open("cipher.txt");
if(decfile.is_open()){

   if (decfile.good()){
       getline(decfile, cipha);

    }
}
decfile.close();

char cstri[cipha.size() + 1];
	strcpy(cstri, cipha.c_str());

 //Initialize the number decry
  mpz_init(decry);
  mpz_set_ui(decry,0);

//Parse the input string as a base 10 number
  mpz_set_str(decry,cstri, 10);
  mpz_class resul (decry);

  ifstream infile;
string line1,line2,line3,line4,line5,line6;
mpz_t num1,num2,num3,num4,num5,num6;


//reading all parameters from the text file
infile.open("parameters.txt");
if(infile.is_open()){

   if (infile.good()){
       getline(infile, line1);
       getline(infile, line2);
       getline(infile, line3);
       getline(infile, line4);
       getline(infile, line5);
       getline(infile, line6);
    }
}
infile.close();
//captures p
char mestr1[line1.size() + 1];
strcpy(mestr1, line1.c_str());
  mpz_init(num1);
  mpz_set_ui(num1,0);
  mpz_set_str(num1,mestr1, 10);
  mpz_class a (num1);
//captures q
  char mestr2[line2.size() + 1];
strcpy(mestr2, line2.c_str());
  mpz_init(num2);
  mpz_set_ui(num2,0);
  mpz_set_str(num2,mestr2, 10);
  mpz_class b(num2);
//captures phi
  char mestr3[line3.size() + 1];
strcpy(mestr3, line3.c_str());
  mpz_init(num3);
  mpz_set_ui(num3,0);
  mpz_set_str(num3,mestr3, 10);
  mpz_class c (num3);
//captures n
  char mestr4[line4.size() + 1];
strcpy(mestr4, line4.c_str());
  mpz_init(num4);
  mpz_set_ui(num4,0);
  mpz_set_str(num4,mestr4, 10);
  mpz_class d (num4);
//captures e
  char mestr5[line5.size() + 1];
strcpy(mestr5, line5.c_str());
  mpz_init(num5);
  mpz_set_ui(num5,0);
  mpz_set_str(num5,mestr5, 10);
  mpz_class e (num5);
//captures d
  char mestr6[line6.size() + 1];
strcpy(mestr6, line6.c_str());
  mpz_init(num6);
  mpz_set_ui(num6,0);
  mpz_set_str(num6,mestr6, 10);
  mpz_class f (num6);

//we shall calculate the message
mpz_t messag;
mpz_init(messag);
mpz_set_ui(messag,0);
mpz_class  mess (messag);
//m=c power d mod n
mpz_powm (mess.get_mpz_t(), resul.get_mpz_t(), f.get_mpz_t(), d.get_mpz_t());
//store results in message.txt
  ofstream messoutfile;
  messoutfile.open("message.txt");
  messoutfile << mess.get_mpz_t()<<endl ;
  messoutfile.close();
  cout<<"Decryption complete and stored in message.txt"<<endl;
  return 0;
}
}
else//error if user did not select 1,2 or 3
{
    printf ("Invalid selection. ");
    return 0;
}

}
