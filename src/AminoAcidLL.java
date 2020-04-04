class AminoAcidLL{
  private char aminoAcid;
  private String[] codons;
  private int[] counts;
  private AminoAcidLL next;

  AminoAcidLL(){
  }

  /********************************************************************************************/
  /* Creates a new node, with a given amino acid/codon 
   * pair and increments the codon counter for that codon.
   * NOTE: Does not check for  repeats!! */
  AminoAcidLL(String inCodon){
      aminoAcid = AminoAcidResources.getAminoAcidFromCodon(inCodon);
      codons = AminoAcidResources.getCodonListForAminoAcid(aminoAcid);
      counts = new int[codons.length];
      addCodon(inCodon);
  }

  /********************************************************************************************/
  /* Recursive method that increments the count for a specific codon:
   * If it should be at this node, increments it and stops, 
   * if not passes the task to the next node. 
   * If there is no next node, add a new node to the list that would contain the codon. 
   */
  private void addCodon(String inCodon){
      if (AminoAcidResources.getAminoAcidFromCodon(inCodon) == aminoAcid) {
        incrementCodon(inCodon);  //Use helper method
      }
      else if (next != null)
        next.addCodon(inCodon);  //Check node
      else {
        next = new AminoAcidLL(inCodon);
      }
  }
  //HELPER METHOD ADDED BY ME
  private void incrementCodon(String inCodon){
    for (int i = 0; i < codons.length; i++) {
      if ( inCodon.equals(codons[i]) )
          counts[i]++;  //Change counts and not codons
    }
  }


  /********************************************************************************************/
  /* Shortcut to find the total number of instances of this amino acid */
  private int totalCount(){
    int sum = 0;

    for (int i = 0; i < counts.length; i++)
      sum += counts[i];

    return sum;
  }

  /********************************************************************************************/
  /* helper method for finding the list difference on two matching nodes
  *  must be matching, but this is not tracked */
  private int totalDiff(AminoAcidLL inList){
    return Math.abs(totalCount() - inList.totalCount());
  }


  /********************************************************************************************/
  /* helper method for finding the list difference on two matching nodes
  *  must be matching, but this is not tracked */
  //TO BE USED
  private int codonDiff(AminoAcidLL inList){
    int diff = 0;
    for(int i=0; i<codons.length; i++){
      diff += Math.abs(counts[i] - inList.counts[i]);
    }
    return diff;
  }

  /********************************************************************************************/
  /* Recursive method that finds the differences in **Amino Acid** counts. 
   * the list *must* be sorted to use this method */
  public int aminoAcidCompare(AminoAcidLL inList){
        if (inList == null && this == null) return 0;

        if (inList == null) {
            if (next == null) return totalCount();
            return totalCount() + next.aminoAcidCompare(null);
        }

        if (next == null) {
            if (inList.next == null)  return totalDiff(inList);
            int a = totalCount();
            while (inList != null) {
               a -= inList.totalCount();
               inList = inList.next;
            }
            return a;
        }

        if (aminoAcid == inList.aminoAcid)
            return totalDiff(inList) + next.aminoAcidCompare(inList.next);

        if (aminoAcid < inList.aminoAcid)
            return totalCount() + next.aminoAcidCompare(inList);

        if (aminoAcid > inList.aminoAcid)
            return - inList.totalCount() + next.aminoAcidCompare(inList);

        return -1;
  }

  /********************************************************************************************/
  /* Same as above, but counts the codon usage differences
   * Must be sorted. */
  public int codonCompare(AminoAcidLL inList){
      if (inList == null && this == null) return 0;

      if (inList == null) {
          if (next == null) return totalCount();
          return totalCount() + next.codonCompare(null);
      }

      if (next == null) {
          if (inList.next == null) return codonDiff(inList);
          int a = totalCount();
          while (inList != null) {
              a -= inList.totalCount();
              inList = inList.next;
          }
          return a;
      }

      if (aminoAcid == inList.aminoAcid)
          return codonDiff(inList) + next.codonCompare(inList.next);

      if (aminoAcid < inList.aminoAcid)
          return totalCount() + next.codonCompare(inList);

      if (aminoAcid > inList.aminoAcid)
          return - inList.totalCount() + next.codonCompare(inList);

      return -1;

  }

  /********************************************************************************************/
  /* Recursively returns the total list of amino acids in the order
   * that they are in in the linked list. */
  public char[] aminoAcidList(){
      if (next == null) {
          char[] arr = {aminoAcid};
          return arr;
      }
      char[] arr = next.aminoAcidList();
      char[] ret = new char[arr.length + 1];
      ret[0] = aminoAcid;
      for (int i = 0; i < arr.length; i++)
          ret[i+1] = arr[i];
      return ret;
  }

  /********************************************************************************************/
  /* Recursively returns the total counts of amino acids in the order that they are in in the linked list. */
  public int[] aminoAcidCounts(){
    if (next == null) {
        int[] arr = {totalCount()};
        return arr;
    }
    int[] arr = next.aminoAcidCounts();
    int[] ret = new int[arr.length + 1];
    ret[0] = totalCount();
    for (int i = 0; i < arr.length; i++)
        ret[i+1] = arr[i];
    return ret;
  }

  /********************************************************************************************/
  /* recursively determines if a linked list is sorted or not */
  public boolean isSorted(){
    if (this.next == null)  return true;
    if (aminoAcid < next.aminoAcid) return next.isSorted();
    return false;
  }

  /********************************************************************************************/
  /* Static method for generating a linked list from an RNA sequence */
  public static AminoAcidLL createFromRNASequence(String inSequence){
    if (inSequence.length() < 3) return null;
    else if (inSequence.length() < 6) return new AminoAcidLL(inSequence.substring(0,3));
    else {
        AminoAcidLL ret = new AminoAcidLL(inSequence.substring(0,3));
        int i = 3;
        while (i + 3 <= inSequence.length()) {
            ret.addCodon(inSequence.substring(i, i + 3));
            i += 3;
        }
        return ret;
    }
  }

  /********************************************************************************************/
  /* sorts a list by amino acid character*/
  public static AminoAcidLL sort(AminoAcidLL inList){
    if (inList == null || inList.next == null)
        return inList;

    AminoAcidLL middle = getMiddle(inList);
    AminoAcidLL middleNext = middle.next;

    middle.next = null;

    AminoAcidLL left = sort(inList);

    AminoAcidLL right = sort(middleNext);

    AminoAcidLL sortedList = sortedMerge(left, right);
    return sortedList;
  }

  public static AminoAcidLL getMiddle(AminoAcidLL head) {
    if (head == null) return null;

    AminoAcidLL slow = head;
    AminoAcidLL fast = head;

    while ( fast.next != null && fast.next.next != null) {
        slow = slow.next;
        fast = fast.next.next;
    }

    return slow;
  }

  public static AminoAcidLL sortedMerge(AminoAcidLL left, AminoAcidLL right) {
    AminoAcidLL sortedList = null;

    if ( left == null ) return right;
    if ( right == null ) return left;

    if (left.aminoAcid <= right.aminoAcid) {
        sortedList = left;
        sortedList.next = sortedMerge(left.next, right);
    }
    else {
        sortedList = right;
        sortedList.next = sortedMerge(left, right.next);
    }

    return sortedList;
  }

  /********************************************************************************************/
  /* Checks if an AminoAcid Linked List is equal to other */
  public boolean AALLCompare(AminoAcidLL inList) {
      if (this == null && inList == null) return true;
      if (next == null && inList.next == null) {
          if (singleCodonCompare(inList)) return true;
          else return false;
      }
      if (!singleCodonCompare(inList)) return false;
      return next.AALLCompare(inList.next);
  }

  public boolean singleCodonCompare (AminoAcidLL inList) {
      if (aminoAcid != inList.aminoAcid) return false;
      for (int i = 0; i < counts.length; i++)
          if (counts[i] != inList.counts[i]) return false;
      return true;
  }

  /*Prints AminoAcid Linked List Character and codons with its counts*/
  public static void print(AminoAcidLL list) {
    while (list != null) {
       System.out.println("\n" + list.aminoAcid);
       for (int i = 0; i < list.codons.length; i++) {
          System.out.print(list.codons[i] + "\t");
       }
       System.out.println();
       for (int i = 0; i < list.counts.length; i++) {
          System.out.print(list.counts[i] + "\t");
       }
       list = list.next;
    }
  }

  /*
  //ELIMINATE.JUST FOR TESTING
  public static void print1(char[] a) {
      for (int i = 0; i < a.length; i++)
          System.out.print(a[i] + " ");
  }
  //ELIMINATE.JUST FOR TESTING
  public static void print2(int[] a) {
      for (int i = 0; i < a.length; i++)
          System.out.print(a[i] + " ");
  } */
}