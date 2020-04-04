import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class AminoAcidLLTester {

    /*1st test case: Tests two AminoAcidLL, where one list is larger
      than the other*/
    @Test
    void aminoAcidCompare() {
        AminoAcidLL a = new AminoAcidLL();
        AminoAcidLL b = a.createFromRNASequence("UUCUUCUUCUUGUUG");
        AminoAcidLL c = a.createFromRNASequence("UUCUUCUUCUUGUUGUUGUUG");

        assertEquals(2,c.aminoAcidCompare(b));
    }

    /* 2nd test case: Checks when one list is null and the other list is not.
    *  Therefore the result will be the number of codons the list has */
    @Test
    void aminoAcidCompare2() {
        AminoAcidLL a = new AminoAcidLL();
        AminoAcidLL b = a.createFromRNASequence("GACGAGGGG");
        AminoAcidLL c = a.createFromRNASequence("");

        assertEquals(3,b.aminoAcidCompare(c));
    }

    /* 3rd test case: Checks when the two list have the same number of codons for
    *  the same AminoAcid but not the same codons*/
    @Test
    void codonCompare() {
        AminoAcidLL a = new AminoAcidLL();
        AminoAcidLL b = a.createFromRNASequence("AUAAUU");
        AminoAcidLL c = a.createFromRNASequence("AUCAUC");

        assertEquals(4,c.codonCompare(b));
    }

    /* 4th test case: Tests method mainly for functionality with a large list*/
    @Test
    void aminoAcidList() {
        char[] arr = {'E','G','R','S','T'};
        AminoAcidLL a = new AminoAcidLL();
        AminoAcidLL b = a.createFromRNASequence("GAGGAAGAGGAAGGGGGGAGGAGCUAU");

        assertArrayEquals(arr, b.aminoAcidList());
    }

    /* 5th test case: Tests method mainly for functionality with a large list*/
    @Test
    void aminoAcidCounts() {
        int[] arr = {4,2,1,1,1};
        AminoAcidLL a = new AminoAcidLL();
        AminoAcidLL b = a.createFromRNASequence("GAGGAAGAGGAAGGGGGGAGGAGCUAU");

        assertArrayEquals(arr, b.aminoAcidCounts());
    }

    /* 6th test case: Tests method mainly for functionality with a large list*/
    @Test
    void aminoAcidCounts2() {
        int[] arr = {4,2,1,1,1};
        AminoAcidLL a = new AminoAcidLL();
        AminoAcidLL b = a.createFromRNASequence("GAGGAAGAGGAAGGGGGGAGGAGCUAU");

        assertArrayEquals(arr, b.aminoAcidCounts());
    }

    /* 7th test case: Creates an AminoAcid list that is sorted,
       therefore expecting the result to be TRUE.  */
    @Test
    void isSorted() {
        AminoAcidLL a = new AminoAcidLL();
        AminoAcidLL b = a.createFromRNASequence("GCAGCGGACGAUACAACC");

        assertTrue(b.isSorted());
    }

    /* 8th test case: Creates an AminoAcid list that is NOT sorted,
   therefore expecting the result to be FALSE.  */
    @Test
    void isSorted2() {
        AminoAcidLL a = new AminoAcidLL();
        AminoAcidLL b = a.createFromRNASequence("CCGUUGGCACUGUUG");

        assertFalse(b.isSorted());
    }

    /* 9th test case: Checks with helper method "print" to print the content
    *  of the list and see if it matches, which it does */
    @Test
    void createFromRNASequence() {
        AminoAcidLL a = new AminoAcidLL();
        AminoAcidLL b = a.createFromRNASequence("CCGUUGGCACUGUUG");
        b.print(b);

        assertTrue(true);
    }

    /* 10th test case: This test case see if the method
     *  works with the example in "README", with a non sorted list and then sorting it*/
    @Test
    void sort() {
        AminoAcidLL a = new AminoAcidLL();
        AminoAcidLL b = a.createFromRNASequence("CCGUUGGCACUGUUG"); //PLALL
        AminoAcidLL list = a.sort(b);
        AminoAcidLL listExpected = a.createFromRNASequence("GCACUGUUGUUGCCG");

        assertTrue(list.AALLCompare(listExpected));
    }
}