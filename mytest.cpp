// CMSC 341 - Spring 2025 - Project 4
#include "dnadb.h"
#include <math.h>
#include <algorithm>
#include <random>
#include <vector>
using namespace std;

unsigned int hashCode(const string str);
string sequencer(int size, int seedNum);

enum RANDOM {UNIFORMINT, UNIFORMREAL, NORMAL, SHUFFLE};
class Random {
public:
    Random(){}
    Random(int min, int max, RANDOM type=UNIFORMINT, int mean=50, int stdev=20) : m_min(min), m_max(max), m_type(type)
    {
        if (type == NORMAL){
            //the case of NORMAL to generate integer numbers with normal distribution
            m_generator = std::mt19937(m_device());
            //the data set will have the mean of 50 (default) and standard deviation of 20 (default)
            //the mean and standard deviation can change by passing new values to constructor 
            m_normdist = std::normal_distribution<>(mean,stdev);
        }
        else if (type == UNIFORMINT) {
            //the case of UNIFORMINT to generate integer numbers
            // Using a fixed seed value generates always the same sequence
            // of pseudorandom numbers, e.g. reproducing scientific experiments
            // here it helps us with testing since the same sequence repeats
            m_generator = std::mt19937(10);// 10 is the fixed seed value
            m_unidist = std::uniform_int_distribution<>(min,max);
        }
        else if (type == UNIFORMREAL) { //the case of UNIFORMREAL to generate real numbers
            m_generator = std::mt19937(10);// 10 is the fixed seed value
            m_uniReal = std::uniform_real_distribution<double>((double)min,(double)max);
        }
        else { //the case of SHUFFLE to generate every number only once
            m_generator = std::mt19937(m_device());
        }
    }
    void setSeed(int seedNum){
        // we have set a default value for seed in constructor
        // we can change the seed by calling this function after constructor call
        // this gives us more randomness
        m_generator = std::mt19937(seedNum);
    }
    void init(int min, int max){
        m_min = min;
        m_max = max;
        m_type = UNIFORMINT;
        m_generator = std::mt19937(10);// 10 is the fixed seed value
        m_unidist = std::uniform_int_distribution<>(min,max);
    }
    void getShuffle(vector<int> & array){
        // this function provides a list of all values between min and max
        // in a random order, this function guarantees the uniqueness
        // of every value in the list
        // the user program creates the vector param and passes here
        // here we populate the vector using m_min and m_max
        for (int i = m_min; i<=m_max; i++){
            array.push_back(i);
        }
        shuffle(array.begin(),array.end(),m_generator);
    }

    void getShuffle(int array[]){
        // this function provides a list of all values between min and max
        // in a random order, this function guarantees the uniqueness
        // of every value in the list
        // the param array must be of the size (m_max-m_min+1)
        // the user program creates the array and pass it here
        vector<int> temp;
        for (int i = m_min; i<=m_max; i++){
            temp.push_back(i);
        }
        std::shuffle(temp.begin(), temp.end(), m_generator);
        vector<int>::iterator it;
        int i = 0;
        for (it=temp.begin(); it != temp.end(); it++){
            array[i] = *it;
            i++;
        }
    }

    int getRandNum(){
        // this function returns integer numbers
        // the object must have been initialized to generate integers
        int result = 0;
        if(m_type == NORMAL){
            //returns a random number in a set with normal distribution
            //we limit random numbers by the min and max values
            result = m_min - 1;
            while(result < m_min || result > m_max)
                result = m_normdist(m_generator);
        }
        else if (m_type == UNIFORMINT){
            //this will generate a random number between min and max values
            result = m_unidist(m_generator);
        }
        return result;
    }

    double getRealRandNum(){
        // this function returns real numbers
        // the object must have been initialized to generate real numbers
        double result = m_uniReal(m_generator);
        // a trick to return numbers only with two deciaml points
        // for example if result is 15.0378, function returns 15.03
        // to round up we can use ceil function instead of floor
        result = std::floor(result*100.0)/100.0;
        return result;
    }

    string getRandString(int size){
        // the parameter size specifies the length of string we ask for
        // to use ASCII char the number range in constructor must be set to 97 - 122
        // and the Random type must be UNIFORMINT (it is default in constructor)
        string output = "";
        for (int i=0;i<size;i++){
            output = output + (char)getRandNum();
        }
        return output;
    }
    
    int getMin(){return m_min;}
    int getMax(){return m_max;}
    private:
    int m_min;
    int m_max;
    RANDOM m_type;
    std::random_device m_device;
    std::mt19937 m_generator;
    std::normal_distribution<> m_normdist;//normal distribution
    std::uniform_int_distribution<> m_unidist;//integer uniform distribution
    std::uniform_real_distribution<double> m_uniReal;//real uniform distribution
};

class Tester{
public:
    bool testInsertionNonColliding(){
        DnaDb dnadb(MINPRIME, hashCode, LINEAR);
        bool result = true;
        

        DNA dna1("AAAAC", 100001, true);
        DNA dna2("AAAAG", 100002, true);
        DNA dna3("AAAAT", 100003, true);
        DNA dna4("AAAAA", 100003, true);
        
        int expectedIndex1 = hashCode(dna1.getSequence()) % MINPRIME;
        int expectedIndex2 = hashCode(dna2.getSequence()) % MINPRIME;
        int expectedIndex3 = hashCode(dna3.getSequence()) % MINPRIME;
        int expectedIndex4 = hashCode(dna4.getSequence()) % MINPRIME;
        if(!dnadb.insert(dna1) || !dnadb.insert(dna2) || !dnadb.insert(dna3) || !dnadb.insert(dna4)){
            result = false;
        }
        if(dnadb.m_currentTable[expectedIndex1]->getSequence() != dna1.getSequence() 
            || dnadb.m_currentTable[expectedIndex2]->getSequence() != dna2.getSequence()
            || dnadb.m_currentTable[expectedIndex3]->getSequence() != dna3.getSequence()
            || dnadb.m_currentTable[expectedIndex4]->getSequence() != dna4.getSequence()){
                result = false;
            }
        
        
        if (dnadb.m_currentSize != 4) {
            result = false;
        }
        return result;
    }
    bool testFindNonExistent(){
        DnaDb dnadb(MINPRIME, hashCode, LINEAR);
        bool result = true;
        DNA dna1("AAAAC", 100001, true);
        DNA dna2("AAAAG", 100002, true);
        DNA dna3("AAAAT", 100003, true);
        DNA dna4("AAAAA", 100003, true);
        dnadb.insert(dna1);
        dnadb.insert(dna2);
        dnadb.insert(dna3);
        dnadb.insert(dna4);
        if(dnadb.getDNA("CCCCC", 100004).getUsed()){
            result = false;
        }
        return result;
    }
    bool testFindNonColliding(){
        DnaDb dnadb(MINPRIME, hashCode, LINEAR);
        bool result = true;
        DNA dna1("AAAAC", 100001, true);
        DNA dna2("AAAAG", 100002, true);
        DNA dna3("AAAAT", 100003, true);
        DNA dna4("AAAAA", 100003, true);
        dnadb.insert(dna1);
        dnadb.insert(dna2);
        dnadb.insert(dna3);
        dnadb.insert(dna4);
        if(!dnadb.getDNA(dna1.getSequence(), dna1.getLocId()).getUsed() 
            || !dnadb.getDNA(dna2.getSequence(), dna2.getLocId()).getUsed()
            || !dnadb.getDNA(dna3.getSequence(), dna3.getLocId()).getUsed()
            || !dnadb.getDNA(dna4.getSequence(), dna4.getLocId()).getUsed()){
            result = false;
        }
        return result;
    }
    bool testFindColliding(){
        DnaDb dnadb(MINPRIME, hashCode, LINEAR);
        bool result = true;
        for (int i = 0; i < floor(MINPRIME * 0.4); i++){
            dnadb.insert(DNA("AAAAA", 100001 + i, true));
            if(!dnadb.getDNA("AAAAA", 100001 + i).getUsed()){
                result = false;
            }
        }
        return result;
    }
    bool testRemoveNonColliding(){
        DnaDb dnadb(MINPRIME, hashCode, LINEAR);
        bool result = true;
        DNA dna1("AAAAC", 100001, true);
        DNA dna2("AAAAG", 100002, true);
        DNA dna3("AAAAT", 100003, true);
        DNA dna4("AAAAA", 100003, true);
        dnadb.insert(dna1);
        dnadb.insert(dna2);
        dnadb.insert(dna3);
        dnadb.insert(dna4);
        
        if(!dnadb.remove(dna1)) {
            result = false;
        }
        if(dnadb.getDNA(dna1.getSequence(), dna1.getLocId()).getUsed()) {
            result = false;
        }
        
        if(!dnadb.remove(dna2)) {
            result = false;
        }
        if(dnadb.getDNA(dna2.getSequence(), dna2.getLocId()).getUsed()) {
            result = false;
        }

        if(!dnadb.getDNA(dna3.getSequence(), dna3.getLocId()).getUsed() ||
           !dnadb.getDNA(dna4.getSequence(), dna4.getLocId()).getUsed()) {
            result = false;
        }
        
        return result;
    }
    bool testRemoveColliding(){
        DnaDb dnadb(MINPRIME, hashCode, LINEAR);
        bool result = true;
        for (int i = 0; i < floor(MINPRIME * 0.4); i++){
            dnadb.insert(DNA("AAAAA", 100001 + i, true));
        }
        //delete every other node
        for (int i = 0; i < floor(MINPRIME * 0.4); i += 2){
            if(!dnadb.remove(DNA("AAAAA", 100001 + i, true))){
                result = false;
            }
        }
        //check if they actually gone
        for (int i = 0; i < floor(MINPRIME * 0.4); i += 2){
            if(dnadb.getDNA("AAAAA", 100001 + i).getUsed()){
                result = false;
            }
        }
        //check if remaining are still there
        for (int i = 1; i < floor(MINPRIME * 0.4); i += 2){
            if(!dnadb.getDNA("AAAAA", 100001 + i).getUsed()){
                result = false;
            }
        }
        
        return result;
    }
    bool testRehashLoadFactor(){
        DnaDb dnadb(MINPRIME, hashCode, LINEAR);
        bool result = true;
        
        cout << "\nTesting Load Factor Rehash:" << endl;
        cout << "Initial table size: " << MINPRIME << endl;
        
        // First insert some items and verify they're accessible
        int initialItems = floor(MINPRIME * 0.4); // Insert enough to be below MAXLOAD
        vector<DNA> insertedDNAs;
        
        for (int i = 0; i < initialItems; i++){
            DNA newDNA("AAAAA", 100001 + i, true);
            if(!dnadb.insert(newDNA)){
                cout << "Error: Failed to insert DNA at index " << i << endl;
                result = false;
            }
            insertedDNAs.push_back(newDNA);
        }
        
        cout << "After initial inserts:" << endl;
        cout << "Current size: " << dnadb.m_currentSize << endl;
        cout << "Current load factor: " << dnadb.lambda() << endl;
        cout << "MAXLOAD: " << MAXLOAD << endl;
        
        // Verify initial items are accessible
        for(int i = 0; i < insertedDNAs.size(); i++){
            if(!dnadb.getDNA(insertedDNAs[i].getSequence(), insertedDNAs[i].getLocId()).getUsed()){
                cout << "Error: DNA at index " << i << " is not accessible after initial inserts" << endl;
                result = false;
            }
        }
        
        // Now insert more items to trigger rehash
        int additionalItems = floor(MINPRIME * 0.2); // This should push us over MAXLOAD
        for (int i = 0; i < additionalItems; i++){
            DNA newDNA("AAAAA", 100001 + initialItems + i, true);
            if(!dnadb.insert(newDNA)){
                cout << "Error: Failed to insert DNA at index " << (initialItems + i) << endl;
                result = false;
            }
            insertedDNAs.push_back(newDNA);
        }
        
        cout << "\nAfter additional inserts:" << endl;
        cout << "Current size: " << dnadb.m_currentSize << endl;
        cout << "Current load factor: " << dnadb.lambda() << endl;
        cout << "Old table exists: " << (dnadb.m_oldTable != nullptr) << endl;
        
        // Verify rehash was triggered
        if(dnadb.m_oldTable == nullptr){
            cout << "Error: Rehash was not triggered despite load factor > MAXLOAD" << endl;
            result = false;
        }
        
        // Insert items until rehash completes
        int finalInserts = 0;
        while(dnadb.m_oldTable != nullptr && finalInserts < 4){
            DNA newDNA("AAAAA", 100001 + initialItems + additionalItems + finalInserts, true);
            if(!dnadb.insert(newDNA)){
                cout << "Error: Failed to insert DNA during rehash completion at index " << (initialItems + additionalItems + finalInserts) << endl;
                result = false;
            }
            insertedDNAs.push_back(newDNA);
            finalInserts++;
            cout << "Final insert " << finalInserts << " completed. Old table still exists: " << (dnadb.m_oldTable != nullptr) << endl;
        }
        
        // Verify rehash completed
        if(dnadb.m_oldTable != nullptr){
            cout << "Error: Rehash did not complete after " << finalInserts << " additional inserts" << endl;
            result = false;
        }
        
        // Verify all items are still accessible and in the new table
        for(int i = 0; i < insertedDNAs.size(); i++){
            if(!dnadb.getDNA(insertedDNAs[i].getSequence(), insertedDNAs[i].getLocId()).getUsed()){
                cout << "Error: DNA at index " << i << " is not accessible after rehash" << endl;
                result = false;
            }
        }
        
        return result;
    }
    bool testRehashDeleteRatio(){
        DnaDb dnadb(MINPRIME, hashCode, LINEAR);
        bool result = true;
        
        // Insert enough items to have a good number of items to delete
        int itemsToInsert = floor(MINPRIME * 0.4); // This will give us enough items to work with
        vector<DNA> insertedDNAs;
        
        cout << "\nTesting Delete Ratio Rehash:" << endl;
        cout << "Initial table size: " << MINPRIME << endl;
        
        // Insert items and store them for later verification
        for (int i = 0; i < itemsToInsert; i++){
            DNA newDNA("AAAAA", 100001 + i, true);
            if(!dnadb.insert(newDNA)){
                cout << "Error: Failed to insert DNA at index " << i << endl;
                result = false;
            }
            insertedDNAs.push_back(newDNA);
        }
        
        cout << "Initial size: " << dnadb.m_currentSize << endl;
        cout << "Initial deleted: " << dnadb.m_currNumDeleted << endl;
        
        // Delete enough items to trigger rehash (delete ratio > 0.8)
        for (int i = 0; i < floor(itemsToInsert * 0.83); i++){
            if(!dnadb.remove(insertedDNAs[i])){
                cout << "Error: Failed to remove DNA at index " << i << endl;
                result = false;
            }
        }
        
        cout << "After deletions - Size: " << dnadb.m_currentSize << endl;
        cout << "After deletions - Deleted: " << dnadb.m_currNumDeleted << endl;
        cout << "Delete ratio: " << dnadb.deletedRatio() << endl;
        cout << "MAXDELETED: " << MAXDELETED << endl;
        
        // Verify rehash was triggered by checking delete ratio and old table
        if(dnadb.m_oldTable == nullptr || dnadb.deletedRatio() <= MAXDELETED){
            cout << "Rehash trigger check failed:" << endl;
            cout << "Old table exists: " << (dnadb.m_oldTable != nullptr) << endl;
            cout << "Delete ratio: " << dnadb.deletedRatio() << endl;
            result = false;
        }
        
        // Delete items until rehash completes
        int additionalDeletes = 0;
        while(dnadb.m_oldTable != nullptr && additionalDeletes < 4){
            int index = floor(itemsToInsert * 0.81) + additionalDeletes;
            if(index < itemsToInsert && !dnadb.remove(insertedDNAs[index])){
                cout << "Error: Failed to remove DNA during rehash completion at index " << index << endl;
                result = false;
            }
            additionalDeletes++;
            cout << "Additional delete " << additionalDeletes << " completed. Old table still exists: " << (dnadb.m_oldTable != nullptr) << endl;
        }
        
        // Verify rehash completed
        if(dnadb.m_oldTable != nullptr){
            cout << "Rehash did not complete after " << additionalDeletes << " additional deletes" << endl;
            result = false;
        }
        
        // Verify remaining items are still accessible and in the new table
        for(int i = floor(itemsToInsert * 0.81) + additionalDeletes; i < itemsToInsert; i++){
            if(!dnadb.getDNA(insertedDNAs[i].getSequence(), insertedDNAs[i].getLocId()).getUsed()){
                cout << "Error: DNA at index " << i << " is not accessible after rehash" << endl;
                result = false;
            }
        }
        
        return result;
    }
};

int main(){
    Tester tester;
    cout << endl;
    cout << "Testing insertion non-colliding case: " << (tester.testInsertionNonColliding() ? "Passed" : "Failed") << endl;
    cout << "Testing find non-existent case: " << (tester.testFindNonExistent() ? "Passed" : "Failed") << endl;
    cout << "Testing find non-colliding case: " << (tester.testFindNonColliding() ? "Passed" : "Failed") << endl;
    cout << "Testing find colliding case: " << (tester.testFindColliding() ? "Passed" : "Failed") << endl;
    cout << "Testing remove non-colliding case: " << (tester.testRemoveNonColliding() ? "Passed" : "Failed") << endl;
    cout << "Testing remove colliding case: " << (tester.testRemoveColliding() ? "Passed" : "Failed") << endl;
    cout << "Testing rehash load factor case: " << (tester.testRehashLoadFactor() ? "Passed" : "Failed") << endl;
    cout << "Testing rehash delete ratio case: " << (tester.testRehashDeleteRatio() ? "Passed" : "Failed") << endl;

    return 0;
}

unsigned int hashCode(const string str) {
    unsigned int val = 0 ;
    const unsigned int thirtyThree = 33 ;  // magic number from textbook
    for ( int i = 0 ; i < str.length(); i++)
       val = val * thirtyThree + str[i] ;
    return val ;
 }
 string sequencer(int size, int seedNum){
     //this function returns a random DNA sequence
     // size param specifies the size of string
     string sequence = "";
     Random rndObject(0,3);
     rndObject.setSeed(seedNum);
     for (int i=0;i<size;i++){
         sequence = sequence + ALPHA[rndObject.getRandNum()];
     }
     return sequence;
 }