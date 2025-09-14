package aligner;

import java.util.Arrays;
import shared.Timer;
import structures.IntList;

/**
 * Implements a WaveFront alignment algorithm for global alignment.
 */
public class WaveFrontAlignerViz implements IDAligner {

	/** Main() passes the args and class to Test to avoid redundant code */
	public static <C extends IDAligner> void main(String[] args) throws Exception {
	    StackTraceElement[] stackTrace = Thread.currentThread().getStackTrace();
		@SuppressWarnings("unchecked")
		Class<C> c=(Class<C>)Class.forName(stackTrace[(stackTrace.length<3 ? 1 : 2)].getClassName());
		Test.testAndPrint(c, args);
	}

    public WaveFrontAlignerViz() {}

	@Override
	public final String name() {return "WaveFrontViz";}
    @Override
    public final float align(byte[] a, byte[] b) {return alignStatic(a, b, null);}
    @Override
    public final float align(byte[] a, byte[] b, int[] pos) {return alignStatic(a, b, pos);}
    @Override
    public final float align(byte[] a, byte[] b, int[] pos, int minScore) {return alignStatic(a, b, pos);}
    @Override
    public final float align(byte[] a, byte[] b, int[] pos, int rStart, int rStop) {
        return alignStatic(a, b, pos, rStart, rStop);
    }
    
    @Override
    public long loops() {return loops;}
	public void setLoops(long x) {loops=x;}

    // Operation types
    private static final byte OP_NONE = 0;
    private static final byte OP_MATCH = 1;
    private static final byte OP_SUB = 2;
    private static final byte OP_INS = 3;
    private static final byte OP_DEL = 4;

    /**
     * Global alignment using WaveFront algorithm with detailed traceback.
     */
    public static final float alignStatic(byte[] query, byte[] ref, int[] posVector) {
        boolean swapped = false;
        
        // Swap to ensure query is not longer than ref
        if(posVector==null && query.length>ref.length) {
            byte[] temp=query;
            query=ref;
            ref=temp;
            swapped = true;
        }

        final int qLen=query.length;
        final int rLen=ref.length;
		long mloops=0;
        
        if(DEBUG_MODE) {
            System.out.println("Query: " + new String(query));
            System.out.println("Ref:   " + new String(ref));
        }
        
        // For visualization
        Visualizer viz=(output==null ? null : new Visualizer(output, 32, 0));
        
//        // Reset loop counter
//        loops = 0;
        
        // Special case for empty sequences
        if(qLen == 0) {
            if(posVector != null) {posVector[0]=0; posVector[1]=Math.max(0, rLen-1);}
            return 0.0f;
        }
        
        // Initialize wavefront data structures
        final int numDiagonals = rLen + qLen + 1;
        final int centerDiagonal = qLen;
        
        // Wavefront arrays
        int[] wf_curr = new int[numDiagonals];
        int[] wf_next = new int[numDiagonals];
        
        // Edit distance matrix
        int[][] editDistance = new int[qLen+1][rLen+1];
        
        // Backtracking matrix
        byte[][] operations = new byte[qLen+1][rLen+1];
        
        // Temporary array for visualization
        int[] rowEditDist = new int[rLen+1];
        
        // Initialize edit distances to infinity
        for(int i=0; i<=qLen; i++) {
            Arrays.fill(editDistance[i], Integer.MAX_VALUE);
        }
        
        // Initialize wavefronts to -1 (not reached)
        Arrays.fill(wf_curr, -1);
        Arrays.fill(wf_next, -1);
        
        // Starting point (0,0)
        wf_curr[centerDiagonal] = 0;
        editDistance[0][0] = 0; // No edits at origin
        
        // For wavefront visualization
        IntList activeList = new IntList(numDiagonals);
        
        // Process wavefronts for increasing edit distances
        for(int distance=0; distance<=qLen+rLen; distance++) {
            if(DEBUG_MODE && distance > 0) {
                System.out.println("WaveFront for edit distance " + distance);
            }
            
            activeList.clear();
            
            // Extend matches along all active diagonals
            for(int diagIdx=0; diagIdx<numDiagonals; diagIdx++) {
                final int diag = diagIdx - centerDiagonal;
                int i = wf_curr[diagIdx];
                
                if(i < 0) continue; // Skip inactive diagonals
                
                int j = i + diag;
                
                // Skip invalid positions
                if(j < 0 || j > rLen) continue;
                
                // Track processed cells
                activeList.add(j);
                loops++;
                
                if(DEBUG_MODE) {
                    System.out.println("Processing diagonal " + diag + " (i=" + i + ", j=" + j + ")");
                }
                
                // Extend matches along this diagonal
                int i_ext = i;
                int j_ext = j;
                
                // Extend exact matches as far as possible
                while(i_ext < qLen && j_ext < rLen && query[i_ext] == ref[j_ext]) {
                    // Record matching extension
                    i_ext++;
                    j_ext++;
                    editDistance[i_ext][j_ext] = distance;
                    operations[i_ext][j_ext] = OP_MATCH;
                    loops++;
                }
                
                // Update position after extension
                if(i_ext > i) {
                    wf_curr[diagIdx] = i_ext;
                    
                    if(DEBUG_MODE) {
                        System.out.println("  Extended to (i=" + i_ext + ", j=" + j_ext + ")");
                    }
                }
                
                // Check if we've reached the end
                if(i_ext == qLen && j_ext == rLen) {
                    if(DEBUG_MODE) {
                        System.out.println("Found complete alignment with edit distance " + distance);
                    }
                    break;
                }
                
                // Next wavefront distance
                int nextDistance = distance + 1;
                
                // Generate next wavefront positions
                
                // Substitution (diagonal move) - edit distance + 1
                if(i_ext < qLen && j_ext < rLen) {
                    if(editDistance[i_ext+1][j_ext+1] > nextDistance) {
                        editDistance[i_ext+1][j_ext+1] = nextDistance;
                        operations[i_ext+1][j_ext+1] = OP_SUB;
                        wf_next[diagIdx] = Math.max(wf_next[diagIdx], i_ext+1);
                    }
                }
                
                // Insertion (gap in reference) - edit distance + 1
                if(i_ext < qLen) {
                    if(editDistance[i_ext+1][j_ext] > nextDistance) {
                        editDistance[i_ext+1][j_ext] = nextDistance;
                        operations[i_ext+1][j_ext] = OP_INS;
                        wf_next[diagIdx-1] = Math.max(wf_next[diagIdx-1], i_ext+1);
                    }
                }
                
                // Deletion (gap in query) - edit distance + 1
                if(j_ext < rLen) {
                    if(editDistance[i_ext][j_ext+1] > nextDistance) {
                        editDistance[i_ext][j_ext+1] = nextDistance;
                        operations[i_ext][j_ext+1] = OP_DEL;
                        wf_next[diagIdx+1] = Math.max(wf_next[diagIdx+1], i_ext);
                    }
                }
            }
            
            // Visualize the current state of the edit distance matrix
            if(PRINT_STEPS && viz != null) {
                // For each row of the reference, extract the edit distances
                for(int j=0; j<=rLen; j++) {
                    // Extract edit distances across this reference position
                    for(int i=0; i<=qLen; i++) {
                        int editDist = editDistance[i][j];
                        rowEditDist[i] = editDist;
                    }
                    // Visualize this row
                    viz.printEditDist(rowEditDist, qLen);
                }
                // Add a separator between distance iterations
                if(viz != null) {
                    viz.print(new IntList(0), rLen, -1);
                }
            }
            
            // If we've found the global alignment, we're done
            if(editDistance[qLen][rLen] != Integer.MAX_VALUE) break;
            
            // Swap wavefronts and reset next wavefront
            int[] tmp = wf_curr;
            wf_curr = wf_next;
            wf_next = tmp;
            Arrays.fill(wf_next, -1);
        }
        
        // Set position vector - for global alignment
        if(posVector != null) {
            posVector[0] = 0;
            posVector[1] = rLen - 1;
        }
        
        // Trace back through the operation matrix to count alignment statistics
        int matches = 0;
        int mismatches = 0;
        int insertions = 0;
        int deletions = 0;
        
        // Use the operation matrix to get the alignment
        StringBuilder alignQ = new StringBuilder();
        StringBuilder alignM = new StringBuilder();
        StringBuilder alignR = new StringBuilder();
        
        int i = qLen;
        int j = rLen;
        
        // Traceback from end to beginning
        while(i > 0 || j > 0) {
            byte op = (i > 0 && j > 0) ? operations[i][j] : OP_NONE;
            
            if(op == OP_MATCH) {
                alignQ.insert(0, (char)query[i-1]);
                alignR.insert(0, (char)ref[j-1]);
                alignM.insert(0, '|');
                matches++;
                i--; j--;
            } else if(op == OP_SUB) {
                alignQ.insert(0, (char)query[i-1]);
                alignR.insert(0, (char)ref[j-1]);
                alignM.insert(0, ' ');
                mismatches++;
                i--; j--;
            } else if(op == OP_INS) {
                alignQ.insert(0, (char)query[i-1]);
                alignR.insert(0, '-');
                alignM.insert(0, ' ');
                insertions++;
                i--;
            } else if(op == OP_DEL) {
                alignQ.insert(0, '-');
                alignR.insert(0, (char)ref[j-1]);
                alignM.insert(0, ' ');
                deletions++;
                j--;
            } else {
                // Default case - if operation not set, use best available move
                if(i > 0 && j > 0) {
                    if(query[i-1] == ref[j-1]) {
                        // Match
                        alignQ.insert(0, (char)query[i-1]);
                        alignR.insert(0, (char)ref[j-1]);
                        alignM.insert(0, '|');
                        matches++;
                    } else {
                        // Mismatch
                        alignQ.insert(0, (char)query[i-1]);
                        alignR.insert(0, (char)ref[j-1]);
                        alignM.insert(0, ' ');
                        mismatches++;
                    }
                    i--; j--;
                } else if(i > 0) {
                    // Insertion
                    alignQ.insert(0, (char)query[i-1]);
                    alignR.insert(0, '-');
                    alignM.insert(0, ' ');
                    insertions++;
                    i--;
                } else {
                    // Deletion
                    alignQ.insert(0, '-');
                    alignR.insert(0, (char)ref[j-1]);
                    alignM.insert(0, ' ');
                    deletions++;
                    j--;
                }
            }
        }

        // Store alignment for debugging
        lastAlignment = "Query: " + alignQ.toString() + "\n" + 
                       "       " + alignM.toString() + "\n" + 
                       "Ref:   " + alignR.toString() + "\n" +
                       "Stats: matches=" + matches + ", mismatches=" + mismatches + 
                       ", insertions=" + insertions + ", deletions=" + deletions;
        
        // If sequences were swapped, swap operations back
        if(swapped) {
            int temp = insertions;
            insertions = deletions;
            deletions = temp;
        }
        
        // Clean up visualization
        if(viz != null) {
        	if(!PRINT_STEPS) {
        		for(int[] row : editDistance) {
        			viz.printEditDist(row, qLen);
        		}
        	}
            viz.shutdown();
        }
        
        // Calculate identity 
        float identity = (float)matches/(float)(matches+mismatches+insertions+deletions);
        
        if(DEBUG_MODE) {
            System.out.println(lastAlignment);
            System.out.println("Final identity: " + identity);
        }
        
        return identity;
    }
    
    /**
     * Wrapper for aligning within a reference window.
     */
    public static final float alignStatic(final byte[] query, final byte[] ref, 
            final int[] posVector, int refStart, int refEnd) {
        refStart = Math.max(refStart, 0);
        refEnd = Math.min(refEnd, ref.length-1);
        final int rlen = refEnd - refStart + 1;
        final byte[] region = (rlen == ref.length ? ref : Arrays.copyOfRange(ref, refStart, refEnd+1));
        final float id = alignStatic(query, region, posVector);
        if(posVector != null) {
            posVector[0] += refStart;
            posVector[1] += refStart;
        }
        return id;
    }
    
    // Debug mode - set to true for detailed alignment output
    private static final boolean DEBUG_MODE = false;
    private static final boolean PRINT_STEPS = false;
    
    // Last alignment result for debugging
    private static String lastAlignment = "";
    
    // For tracking cells processed
    static long loops = 0;
    public static String output = null;
}