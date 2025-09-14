/**
 * Visualizer for alignment matrices that generates ASCII representations of
 * alignment exploration patterns and score distributions.
 * <p>
 * Produces text files showing how the Bridge Building Aligner (BBA) explores
 * the dynamic programming matrix, highlighting explored cells, their scores,
 * and the optimal alignment path.
 *
 * @author Brian Bushnell
 * @date April 2024
 */
package aligner;

import fileIO.ByteStreamWriter;
import shared.Tools;
import structures.ByteBuilder;
import structures.IntList;

public class Visualizer {
    
    /**
     * Creates a new visualization output for alignment matrices.
     * 
     * @param fname Output filename for the visualization
     * @param pBits Number of bits used for position information in the score encoding
     * @param cBits Number of bits used for deletion information in the score encoding
     */
    public Visualizer(String fname, int pBits, int cBits) {
        positionBits=pBits;
        countBits=cBits;
        scoreShift=positionBits+countBits;
        bsw=ByteStreamWriter.makeBSW(fname, true, false, true);
    }
    
    /**
     * Properly terminates the output file.
     * Should be called when visualization is complete.
     */
    public void shutdown() {
        bsw.poisonAndWait();
    }
    
    /**
     * Generates a simple visualization line marking explored cells.
     * Original visualization style that represents exploration pattern
     * without score information.
     * 
     * @param active List of positions that were actively explored
     * @param length Length of the current row
     * @param maxPos Position with the highest score in this row (-1 if none)
     */
    public void print(IntList active, int length, int maxPos) {
        ByteBuilder bb=new ByteBuilder(length);
        
        //Prefill line with spaces
        for(int i=0; i<length; i++) {bb.append(' ');}
        
        //Mark explored cells
        for(int i=0; i<active.size; i++) {
            int cell=active.get(i);
            bb.set(cell, '@');
        }
        
        //Mark the best path
        if(maxPos>=0) {bb.set(maxPos, '*');}
        bb.nl();
        bsw.print(bb);
    }
    
    /** 
     * Wrapper for int[].
     */
    public void print(int[] scores, int bandStart, int bandEnd, int rLen) {
    	long[] scores2=new long[scores.length];
    	for(int i=0; i<scores.length; i++) {scores2[i]=scores[i];}
    	print(scores2, bandStart, bandEnd, rLen);
    }
    
    /** 
     * Wrapper for byte[]
     */
	public void print(byte[] scores, int bandStart, int bandEnd, int rLen) {
    	long[] scores2=new long[scores.length];
    	for(int i=0; i<scores.length; i++) {scores2[i]=scores[i];}
    	print(scores2, bandStart, bandEnd, rLen);
	}
    
    /** 
     * Visualizer for banded aligners.
     * @param scores Array of scores for the current row
     * @param bandStart Min scored position for this row
     * @param bandEnd Max scored position for this row
     * @param rLen Reference length
     */
    public void print(long[] scores, int bandStart, int bandEnd, int rLen) {
    	  ByteBuilder bb=new ByteBuilder(scores.length);
          final int maxPos=Tools.maxIndex(scores);
          final int maxScore=(int)(scores[maxPos]>>scoreShift);
          int minScore=maxScore;
          for(int i=bandStart; i<=bandEnd; i++) {
        	  if(scores[i]>BAD/2) {
        		  minScore=(int)Math.min(minScore, scores[i]>>scoreShift);
        	  }
          }
          
          //Fill good cells score symbols and bad cells with spaces
          for(int i=0; i<scores.length && i<=rLen; i++) {
              final long rawScore=scores[i];
              final int score=(int)(rawScore>>scoreShift);
              byte symbol=(rawScore<=BAD || i<bandStart || i>bandEnd) ? 
            		  (byte)' ' : scoreToSymbol(score, maxScore);
              bb.append(symbol);
          }
          //Mark the best path
          bb.set(maxPos, '*');
          bb.nl();
          bsw.print(bb);
    }

	public void printMSA(int[] scores, int qLen, int rLen, int pointsMatch) {
  	  ByteBuilder bb=new ByteBuilder(scores.length);
      final int maxPos=Tools.maxIndex(scores);
      final int maxScore=(int)(scores[maxPos]>>scoreShift)/pointsMatch;
      final int bad=(-qLen/2);
      //Fill good cells score symbols and bad cells with spaces
      for(int i=0; i<scores.length && i<=rLen; i++) {
          final int rawScore=scores[i];
          final int score=(int)((rawScore>>scoreShift)/pointsMatch);
          byte symbol=(score<=bad) ? 
        		  (byte)' ' : scoreToSymbol(score, maxScore);
          bb.append(symbol);
//          bb.append(score).tab();
      }
      //Mark the best path
      bb.set(maxPos, '*');
      bb.nl();
      bsw.print(bb);
	}
    
    /**
     * Generates a detailed visualization line showing score distributions.
     * Enhanced visualization that represents scores with different characters,
     * providing richer information about the scoring landscape.
     * 
     * @param scores Array of scores for the current row
     * @param active List of positions that were actively explored (optional, may be null)
     */
    public void print(long[] scores, IntList active, int rLen) {
        ByteBuilder bb=new ByteBuilder(scores.length);
        final int maxPos=Tools.maxIndex(scores);
        final int maxScore=(int)(scores[maxPos]>>scoreShift);
        int minScore=maxScore-1;

        
        for(int i=3; i<scores.length; i++) {
        	if(scores[i]>BAD/2) {
        		minScore=(int)Math.min(minScore, scores[i]>>scoreShift);
        	}
        }
        
        //Fill good cells score symbols and bad cells with spaces
        for(int i=0; i<scores.length && i<=rLen; i++) {
            final long rawScore=scores[i];
            final int score=(int)(rawScore>>scoreShift);
            byte symbol=(rawScore<=BAD ? (byte)' ' : scoreToSymbol(score, maxScore));
//            byte symbol=(rawScore<=BAD ? (byte)' ' : scoreToSymbol(score, maxScore, minScore));
            bb.append(symbol);
        }
        //Set explored bad cells to '.' indicating their score was erased
        if(active!=null) {
            for(int i=0; i<active.size; i++) {
                int pos=active.get(i);
                if(pos<=rLen &&bb.get(pos)==' ') {bb.set(pos, '.');}
            }
        }
        //Mark the best path
        bb.set(maxPos, '*');
        bb.nl();
        bsw.print(bb);
    }
    
    /**
     * Generates a detailed visualization line showing score distributions.
     * Enhanced visualization that represents scores with different characters,
     * providing richer information about the scoring landscape.
     * 
     * @param editDist Edit distances for the row, or very high if unexplored
     * @param rLen Reference length.
     */
    public void printEditDist(int[] editDist, int rLen) {
//    	System.err.println("Called viz with "+Arrays.toString(editDist)+", "+rLen);
        ByteBuilder bb=new ByteBuilder(editDist.length);
        final int bad=Integer.MAX_VALUE/2;
        final int minPos=Tools.minIndex(editDist);
        final int minEdits=(editDist[minPos]);
        int maxEdits=minEdits;
        for(int e : editDist) {
        	if(e<bad && e>maxEdits) {maxEdits=e;}
        }
        
        //Fill good cells score symbols and bad cells with spaces
        for(int i=0; i<editDist.length /*&& i<=rLen*/; i++) {
            final int rawScore=editDist[i];
            final int score=Tools.max(0, (minEdits-rawScore+symbols.length-1));
            byte symbol=(rawScore>=bad ? (byte)' ' : symbols[score]);
            bb.append(symbol);
        }
        //Mark the best path
        if(minPos>=0) {bb.set(minPos, '*');}
        bb.nl();
        bsw.print(bb);
    }
    
    /**
     * Converts a numeric score to a display character.
     * Uses absolute or relative scoring based on the maximum score value.
     * 
     * @param score The score to convert
     * @param maxScore The maximum score in the current row
     * @return A character representing the score's magnitude
     */
    private static byte scoreToSymbol(int score, int maxScore) {
        final int symbol;
        if(useAbsolute && !useRelative) {//Absolute;
            symbol=score;
        }else if(useScaled) {
            int minScore=Math.min(0, maxScore-(int)(symbols.length*0.40f));
            int range=maxScore-minScore;
            float mult=symbols.length/(float)range;
            symbol=Math.round(mult*(score-minScore));
        }else if(useAbsolute && useRelative && maxScore>=0 && maxScore<=36) {//Absolute;
            symbol=score+26;
        }else {//Relative
            symbol=score-maxScore+symbols.length;
        }
        return symbols[Tools.mid(0, symbol, symbols.length-1)];
    }
    
//    private static byte scoreToSymbol(int score, int maxScore, int minScore) {
//        final int symbol;
//        if(useAbsolute && !useRelative) {//Absolute;
//            symbol=score;
//        }else if(useScaled) {
//            int range=maxScore-minScore;
//            float mult=symbols.length/(float)range;
//            symbol=Math.round(mult*(score-minScore));
//        }else if(useAbsolute && useRelative && maxScore>=0 && maxScore<=36) {//Absolute;
//            symbol=score+16;
//        }else {//Relative
//            symbol=score-maxScore+symbols.length;
//        }
//        return symbols[Tools.mid(0, symbol, symbols.length-1)];
//    }
    
	private static int[] makeSymbolMap(byte[] symbols) {
		int[] map=new int[128];
		for(int i=0; i<symbols.length; i++) {
			map[symbols[i]]=i;
		}
		return map;
	}
    
    private ByteStreamWriter bsw;
    
    // Bit field definitions
    private final int positionBits;
    private final int countBits;
    private final int scoreShift;

    /** Value representing invalid or pruned cells */
    private static final long BAD=(Long.MIN_VALUE/2);
    public static boolean useAbsolute=true;
    public static boolean useRelative=true;
    public static boolean useScaled=false;
    
    /** Character set for score visualization (a-z, 0-9, A-Z) from lowest to highest */
    static final byte[] symbols=new byte[] {
            'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 
            'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z',
            '0', '1', '2', '3', '4', '5', '6', '7', '8', '9',
            'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 
            'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z'
        };
    static final int[] symbolMap=makeSymbolMap(symbols);
}