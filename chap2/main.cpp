#include <stdio.h>
// Tror det er lettest å forstå koden om du leser main-funksjonen først!

// "Vanlig" funksjon. Denne får alltid en kopi av verdien som blir sendt inn
// Selv om du gjør endringer på pNumber-variabelen vil ikke det påvirke myInt i main-funksjonen
int addTwo(int pNumber){
	pNumber = pNumber + 2;
	return pNumber;
}

// Denne tar imot en pointer (*-tegn) til en int, dvs en peker til hvor i minnet på maskinen den er lagret (addressen).
void changeNumber(int* pNumber){
	// *pNumber betyr at vi mener verdien pointeren peker til
	*pNumber = 4;

	// Her endrer vi addressen pointeren peker til (Setter den til NULL (Det store intet)).
	// Etter å ha gjort dette har vi ingen referanse til tallet som ble sendt av main-funksjonen lengre
	pNumber = NULL;

	// FEIL måte å gjøre det på. De vi sier her er at vi tar peker-addressen og legger til 2.
	// Dvs. at vi flytter pekeren 2 "plasser" videre i minnet, som betyr at vi ikke har anelse om hva den peker mot lengre (Farlig!)
	// Det kan være riktig å gjøre dette hvis du vet hva som ligger etter addressen i minnet (f. eks. om du har en array/liste)
	//pNumber = pNumber + 2;

	// Denne vil få programmet til å kræsje fordi pekeren er satt til NULL. (Man kan ikke endre på "ingenting")
	//*pNumber = 7;
}

int main()
{
	// Vi lager en int, og setter den til verdien 6
	int myInt = 6;
	// Print til konsoll
	printf( "myInt: %i\n", myInt );

	// Vi bruker addTwo for å få en ny verdi, 2 høyere enn myInt (8).
	// Ingen endringer blir gjort på myInt fordi det er en kopi som blir sendt til addTwo-funksjonen
	int addedInt = addTwo(myInt);
	printf( "addedInt: %i\n", addedInt);
	printf( "myInt: %i\n", myInt );

	// Her sender vi addressen i minnet (&-tegn) hvor myInt ligger lagret.
	// Nå kan changeNumber-funksjonen gjøre endringer rett på myInt-variabelen
	changeNumber(&myInt);
	printf( "myInt: %i\n", myInt );
}

// &myInt = addressen til myInt-variabelen
// int* pNumber = en peker (addresse) til en int
// *pNumber = verdien pekeren peker til