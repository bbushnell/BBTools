package aligner;

import java.nio.ByteBuffer;
import java.util.Arrays;

import jdk.incubator.vector.ByteVector;
import jdk.incubator.vector.IntVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorSpecies;

/**
 * Java implementation of the diagonal-based SIMD alignment algorithm from ksw2.
 * Uses anti-diagonal traversal with a drifting window for efficient processing.
 * 
 * @author Isla (based on ksw2 by Heng Li)
 * @date May 2025
 */
public class DiagonalAligner implements IDAligner {

	/** Main() passes the args and class to Test to avoid redundant code */
	public static <C extends IDAligner> void main(String[] args) throws Exception {
	    StackTraceElement[] stackTrace = Thread.currentThread().getStackTrace();
		@SuppressWarnings("unchecked")
		Class<C> c=(Class<C>)Class.forName(stackTrace[(stackTrace.length<3 ? 1 : 2)].getClassName());
		Test.testAndPrint(c, args);
	}

    private static final VectorSpecies<Byte> BSPECIES = ByteVector.SPECIES_128;
    private static final VectorSpecies<Integer> ISPECIES = IntVector.SPECIES_128;
    private static final int BWIDTH = BSPECIES.length(); // Should be 16
    
    @Override
    public String name() {
        return "Diagonal";
    }

    @Override
    public float align(byte[] query, byte[] target) {
        return align(query, target, null);
    }

    @Override
    public float align(byte[] query, byte[] target, int[] posVector) {
        return align(query, target, posVector, 0);
    }

    @Override
    public float align(byte[] query, byte[] target, int[] posVector, int minScore) {
        return align(query, target, posVector, 0, target.length - 1);
    }

    @Override
    public float align(byte[] query, byte[] target, int[] posVector, int refStart, int refEnd) {
        // Initialize parameters
        int qlen = query.length;
        int tlen = target.length;
        
        // Score parameters (can be made configurable)
        byte matchScore = 2;
        byte mismatchScore = -2;
        byte gapOpen = 4;
        byte gapExtend = 1;
        int w = -1; // Bandwidth: -1 means adaptive
        int zdrop = 200; // Z-drop threshold
        int endBonus = 0;
        boolean scoreOnly = posVector == null;
        
        // Create alignment result object
        ExtzResult ez = new ExtzResult();
        
        // Call the actual alignment function
        diagonalAlign(qlen, query, tlen, target, matchScore, mismatchScore, gapOpen, gapExtend, 
                w, zdrop, endBonus, scoreOnly, ez);
        
        // Process result
        if (posVector != null) {
            posVector[0] = Math.max(0, ez.max_t - ez.max_q);
            posVector[1] = Math.min(tlen-1, ez.max_t);
        }
        
        // Calculate identity
        float matches = ez.score / (float)(matchScore);
        float identity = matches / Math.max(qlen, tlen);
        return identity;
    }
    
    // The main alignment function
    private void diagonalAlign(int qlen, byte[] query, int tlen, byte[] target, 
                              byte matchScore, byte mismatchScore, byte gapOpen, byte gapExtend,
                              int w, int zdrop, int endBonus, boolean scoreOnly, ExtzResult ez) {
        
        if (qlen <= 0 || tlen <= 0) return;
        
        // Adjust bandwidth if needed
        if (w < 0) w = Math.max(tlen, qlen);
        int wl = w, wr = w;
        
        // Calculate sizes and allocate memory
        int tlen_ = (tlen + 15) / 16;
        int n_col_ = Math.min(qlen, tlen);
        n_col_ = ((Math.min(n_col_, w + 1) + 15) / 16) + 1;
        int qlen_ = (qlen + 15) / 16;
        
        // Score limits
        byte max_sc = matchScore;
        byte min_sc = mismatchScore;
        
        // Early return if penalties are too high
        if (-min_sc > 2 * (gapOpen + gapExtend)) return;
        
        // Prepare vectors with constants
        ByteVector zero = ByteVector.zero(BSPECIES);
        ByteVector qVec = ByteVector.broadcast(BSPECIES, gapOpen);
        ByteVector qe2Vec = ByteVector.broadcast(BSPECIES, (byte)((gapOpen + gapExtend) * 2));
        ByteVector flag1Vec = ByteVector.broadcast(BSPECIES, (byte)1);
        ByteVector flag2Vec = ByteVector.broadcast(BSPECIES, (byte)2);
        ByteVector flag8Vec = ByteVector.broadcast(BSPECIES, (byte)0x08);
        ByteVector flag16Vec = ByteVector.broadcast(BSPECIES, (byte)0x10);
        ByteVector matchVec = ByteVector.broadcast(BSPECIES, matchScore);
        ByteVector mismatchVec = ByteVector.broadcast(BSPECIES, mismatchScore);
        ByteVector wildcardVec = ByteVector.broadcast(BSPECIES, (byte)-gapExtend); // N score
        ByteVector maxScoreVec = ByteVector.broadcast(BSPECIES, (byte)(matchScore + (gapOpen + gapExtend) * 2));
        
        // Allocate arrays
        byte[][] arrays = new byte[6][tlen_ * 16];
        byte[] u8 = arrays[0];
        byte[] v8 = arrays[1];
        byte[] x8 = arrays[2];
        byte[] y8 = arrays[3];
        byte[] s8 = arrays[4];
        byte[] sf = arrays[5];
        
        // Prepare query and target
        byte[] qr = new byte[qlen];
        for (int t = 0; t < qlen; ++t) {
            qr[t] = query[qlen - 1 - t];
        }
        System.arraycopy(target, 0, sf, 0, tlen);
        
        // Score tracking
        int[] H = new int[tlen_ * 16];
        Arrays.fill(H, Integer.MIN_VALUE / 2);
        
        // Prepare for traceback
        byte[][] p = null;
        int[] off = null, off_end = null;
        if (!scoreOnly) {
            p = new byte[(qlen + tlen - 1)][n_col_ * 16];
            off = new int[qlen + tlen - 1];
            off_end = new int[qlen + tlen - 1];
        }
        
        // Main diagonal loop
        for (int r = 0, last_st = -1, last_en = -1; r < qlen + tlen - 1; ++r) {
            // Find the boundaries
            int st = 0, en = tlen - 1;
            if (st < r - qlen + 1) st = r - qlen + 1;
            if (en > r) en = r;
            if (st < (r - wr + 1) >> 1) st = (r - wr + 1) >> 1; // ceil
            if (en > (r + wl) >> 1) en = (r + wl) >> 1; // floor
            
            if (st > en) {
                ez.zdropped = true;
                break;
            }
            
            int st0 = st, en0 = en;
            st = (st / 16) * 16;
            en = ((en + 16) / 16) * 16 - 1;
            
            // Set boundary conditions
            byte x1 = 0, v1 = 0;
            if (st > 0) {
                if (st - 1 >= last_st && st - 1 <= last_en) {
                    x1 = x8[st - 1];
                    v1 = v8[st - 1];
                }
            } else {
                x1 = 0;
                v1 = (byte)(r != 0 ? gapOpen : 0);
            }
            
            if (en >= r) {
                y8[r] = 0;
                u8[r] = (byte)(r != 0 ? gapOpen : 0);
            }
            
            // Set scores for the current diagonal
            for (int t = st0; t <= en0; t += 16) {
                int blockEnd = Math.min(t + 16, en0 + 1);
                for (int i = t; i < blockEnd; i++) {
                    byte sq = sf[i];
                    byte st_ = qr[i];
                    boolean mask = (sq == 4 || st_ == 4); // Assuming 4 is wildcard/N
                    boolean match = (sq == st_);
                    s8[i] = mask ? wildcardVec.lane(0) : 
                           (match ? matchVec.lane(0) : mismatchVec.lane(0));
                }
            }
            
            // Process the band with vectors
            ByteVector x1Vec = ByteVector.broadcast(BSPECIES, x1);
            ByteVector v1Vec = ByteVector.broadcast(BSPECIES, v1);
            
            int st_ = st / 16, en_ = en / 16;
            
            if (scoreOnly) {
                // Score-only loop (without traceback info)
                for (int t = st_; t <= en_; ++t) {
                    int offset = t * 16;
                    
                    ByteVector sVec = ByteVector.fromArray(BSPECIES, s8, offset);
                    ByteVector zVec = sVec.add(qe2Vec);
                    
                    ByteVector xt1Vec = ByteVector.fromArray(BSPECIES, x8, offset);
                    byte tmp = offset + 15 < x8.length ? x8[offset + 15] : 0;
                    xt1Vec = xt1Vec.lanewise(VectorOperators.LSHL, 1);
                    xt1Vec = xt1Vec.withLane(0, x1);
                    x1 = tmp;
                    
                    ByteVector vt1Vec = ByteVector.fromArray(BSPECIES, v8, offset);
                    tmp = offset + 15 < v8.length ? v8[offset + 15] : 0;
                    vt1Vec = vt1Vec.lanewise(VectorOperators.LSHL, 1);
                    vt1Vec = vt1Vec.withLane(0, v1);
                    v1 = tmp;
                    
                    ByteVector aVec = xt1Vec.add(vt1Vec);
                    ByteVector utVec = ByteVector.fromArray(BSPECIES, u8, offset);
                    ByteVector bVec = ByteVector.fromArray(BSPECIES, y8, offset).add(utVec);
                    
                    zVec = zVec.max(bVec);
                    zVec = zVec.min(maxScoreVec);
                    
                    ByteVector newUVec = zVec.sub(vt1Vec);
                    newUVec.intoArray(u8, offset);
                    
                    ByteVector newVVec = zVec.sub(utVec);
                    newVVec.intoArray(v8, offset);
                    
                    zVec = zVec.sub(qVec);
                    aVec = aVec.sub(zVec);
                    bVec = bVec.sub(zVec);
                    
                    aVec.max(zero).intoArray(x8, offset);
                    bVec.max(zero).intoArray(y8, offset);
                }
            } else {
                // Loop with traceback info generation
                byte[] pr = p[r];
                off[r] = st;
                off_end[r] = en;
                
                for (int t = st_; t <= en_; ++t) {
                    int offset = t * 16;
                    int prOffset = offset - st;
                    
                    ByteVector sVec = ByteVector.fromArray(BSPECIES, s8, offset);
                    ByteVector zVec = sVec.add(qe2Vec);
                    
                    ByteVector xt1Vec = ByteVector.fromArray(BSPECIES, x8, offset);
                    byte tmp = offset + 15 < x8.length ? x8[offset + 15] : 0;
                    xt1Vec = xt1Vec.lanewise(VectorOperators.LSHL, 1);
                    xt1Vec = xt1Vec.withLane(0, x1);
                    x1 = tmp;
                    
                    ByteVector vt1Vec = ByteVector.fromArray(BSPECIES, v8, offset);
                    tmp = offset + 15 < v8.length ? v8[offset + 15] : 0;
                    vt1Vec = vt1Vec.lanewise(VectorOperators.LSHL, 1);
                    vt1Vec = vt1Vec.withLane(0, v1);
                    v1 = tmp;
                    
                    ByteVector aVec = xt1Vec.add(vt1Vec);
                    ByteVector utVec = ByteVector.fromArray(BSPECIES, u8, offset);
                    ByteVector bVec = ByteVector.fromArray(BSPECIES, y8, offset).add(utVec);
                    
                    // Direction tracking for traceback
                    VectorMask<Byte> dirMask = aVec.compare(VectorOperators.GT, zVec);
                    ByteVector dVec = ByteVector.broadcast(BSPECIES, (byte)0);
                    dVec = dVec.blend(flag1Vec, dirMask);
                    
                    zVec = zVec.max(aVec);
                    dirMask = bVec.compare(VectorOperators.GT, zVec);
                    dVec = dVec.blend(flag2Vec, dirMask);
                    
                    zVec = zVec.max(bVec);
                    zVec = zVec.min(maxScoreVec);
                    
                    ByteVector newUVec = zVec.sub(vt1Vec);
                    newUVec.intoArray(u8, offset);
                    
                    ByteVector newVVec = zVec.sub(utVec);
                    newVVec.intoArray(v8, offset);
                    
                    zVec = zVec.sub(qVec);
                    aVec = aVec.sub(zVec);
                    bVec = bVec.sub(zVec);
                    
                    // Track non-zero values
                    VectorMask<Byte> aMask = aVec.compare(VectorOperators.GT, zero);
                    ByteVector newXVec = aVec.max(zero);
                    newXVec.intoArray(x8, offset);
                    dVec = dVec.blend(flag8Vec, aMask);
                    
                    VectorMask<Byte> bMask = bVec.compare(VectorOperators.GT, zero);
                    ByteVector newYVec = bVec.max(zero);
                    newYVec.intoArray(y8, offset);
                    dVec = dVec.blend(flag16Vec, bMask);
                    
                    // Store traceback info
                    dVec.intoArray(pr, prOffset);
                }
            }
            
            // Find max score
            int maxH = Integer.MIN_VALUE / 2;
            int maxT = en0;
            
            if (r > 0) {
                // Update scores and find maximum
                for (int t = st0; t <= en0; ++t) {
                    if (t > 0) {
                        H[t] = H[t-1] + u8[t] - (gapOpen + gapExtend);
                    } else {
                        H[t] += v8[t] - (gapOpen + gapExtend);
                    }
                    
                    if (H[t] > maxH) {
                        maxH = H[t];
                        maxT = t;
                    }
                }
            } else {
                H[0] = v8[0] - (gapOpen + gapExtend) * 2;
                maxH = H[0];
                maxT = 0;
            }
            
            // Update alignment endpoints
            if (en0 == tlen - 1 && H[en0] > ez.mte) {
                ez.mte = H[en0];
                ez.mte_q = r - en;
            }
            
            if (r - st0 == qlen - 1 && H[st0] > ez.mqe) {
                ez.mqe = H[st0];
                ez.mqe_t = st0;
            }
            
            // Apply Z-drop heuristic
            if (applyZDrop(ez, maxH, r, maxT, zdrop, gapExtend)) {
                break;
            }
            
            if (r == qlen + tlen - 2 && en0 == tlen - 1) {
                ez.score = H[tlen - 1];
            }
            
            last_st = st;
            last_en = en;
        }
        
        // Set alignment endpoints based on max score
        if (ez.score <= 0) {
            ez.max_q = qlen - 1;
            ez.max_t = tlen - 1;
        }
    }
    
    // Z-drop heuristic
    private boolean applyZDrop(ExtzResult ez, int curH, int r, int t, int zdrop, int e) {
        if (curH < ez.max) return false;
        
        if (curH > ez.max) {
            ez.max = curH;
            ez.max_t = t;
            ez.max_q = r - t;
        } else if (t >= ez.max_t && r - t >= ez.max_q) {
            ez.max_t = t;
            ez.max_q = r - t;
        }
        
        if (r - ez.max_q - ez.max_t > zdrop / e) {
            ez.zdropped = true;
            ez.max_drop_r = r;
            ez.max_drop_q = ez.max_q;
            ez.max_drop_t = ez.max_t;
            return true;
        }
        
        return false;
    }
    
    // Result class
    static class ExtzResult {
        int max = 0;
        int mqe = 0, mqe_t = -1;
        int mte = 0, mte_q = -1;
        int score = 0;
        int max_t = -1, max_q = -1;
        int max_drop_r = -1, max_drop_q = -1, max_drop_t = -1;
        boolean zdropped = false;
        boolean reach_end = false;
    }

	@Override
	public long loops() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public void setLoops(long i) {
		// TODO Auto-generated method stub
		
	}
}
