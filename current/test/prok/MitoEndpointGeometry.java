package test.prok;


import dna.AminoAcid;
import parse.LineParser1;
import parse.Parse;

/** Exact source-segment reconstruction and sense-oriented genomic flanks.
 * Segments are the extractor's one-based inclusive transcript-order ledger;
 * windows remain genomic (introns are not silently spliced away).
 * Preserved mitochondrial endpoint experiment: this class supports construction,
 * inference or validation of the experimental endpoint method. The method was not
 * promoted, and the mitochondrial-SSU direction was set aside on October 24 because
 * these genes did not fit the prokaryotic calling approach. Kept for a future
 * mitochondrial calling investigation; ordinary CallGenes never enables it.
 * @author Brian Bushnell, Raiden
 */
public final class MitoEndpointGeometry {

	public MitoEndpointGeometry(String segments, String strand, boolean crosses){
		require(segments!=null && !segments.isEmpty(), "A member needs the original segment ledger");
		final byte[] bytes=segments.getBytes(java.nio.charset.StandardCharsets.US_ASCII);
		final LineParser1 p=new LineParser1(',');p.set(bytes);
		lo=new int[p.terms()];hi=new int[lo.length];reverse=new boolean[lo.length];
		for(int i=0; i<lo.length; i++){
			final String s=p.parseString(i);final int dot=s.indexOf(".."),colon=s.indexOf(':');
			require(dot>0 && colon>dot+2 && colon==s.length()-2 && (s.endsWith(":+") || s.endsWith(":-")), "Invalid recorded segment: "+s);
			lo[i]=Parse.parseInt(s, 0, dot);hi[i]=Parse.parseInt(s, dot+2, colon);reverse[i]=s.charAt(colon+1)=='-';
			require(lo[i]>0 && hi[i]>=lo[i], "Recorded segments must be one-based inclusive: "+s);
		}
		boolean mixed=false;for(boolean r:reverse){mixed|=r!=reverse[0];}
		require(strand.equals(mixed?"MIXED":reverse[0]?"-":"+"), "Recorded strand must agree with every source segment");
		mixedStrand=mixed;crossesOrigin=crosses;
	}

	/** Verify exact spliced gene bases without allocating another gene sequence. */
	public void verify(byte[] genome, byte[] exact, String id){
		require(genome!=null && exact!=null && exact.length>0, "Source parity needs nonempty sequences: "+id);
		int pos=0;
		for(int j=0; j<lo.length; j++){
			require(hi[j]<=genome.length, "Source segment exceeds its contig: "+id);
			for(int i=0, n=hi[j]-lo[j]+1; i<n; i++){
				final byte base=reverse[j]?complement(genome[hi[j]-1-i]):genome[lo[j]-1+i];
				require(pos<exact.length && base==exact[pos], "Genomic reconstruction differs from archived member "+id+" at transcript base "+pos);pos++;
			}
		}
		require(pos==exact.length, "Segment lengths do not conserve archived member: "+id);
	}

	/** Return null for geometry lacking a single ordered genomic span, never a fabricated flank. */
	public Window window(byte[] genome, int flank, int maxSpan){
		require(genome!=null && genome.length>0 && flank>=0 && maxSpan>0, "Genomic window needs bounded source coordinates");
		if(mixedStrand){return null;}
		final int len=genome.length;int start=-1,end=-1,turn=0;
		for(int i=0; i<lo.length; i++){
			require(hi[i]<=len, "Segment exceeds the source contig");
			int a=reverse[i]?len-hi[i]:lo[i]-1,b=reverse[i]?len-lo[i]:hi[i]-1;
			if(i>0 && a+turn<=end){if(!crossesOrigin || turn!=0){return null;}turn=len;}
			a+=turn;b+=turn;if(i==0){start=a;}end=b;
		}
		if((turn!=0)!=crossesOrigin || end-start+1>len || end-start+1>maxSpan){return null;}
		// Only explicitly origin-crossing annotations justify wrapping; ordinary
		// contig ends remain clipped because topology is absent from this ledger.
		final int span=end-start+1,free=len-span;
		final int left=crossesOrigin?Math.min(flank, free/2):Math.min(flank, start);
		final int right=crossesOrigin?Math.min(flank, free-left):Math.min(flank, len-1-end);
		final byte[] bases=new byte[left+span+right];int q=start-left;
		for(int i=0; i<bases.length; i++,q++){
			int at=q;if(at<0){at+=len;}else if(at>=len){at-=len;}
			assert(at>=0 && at<len):"One-turn source windows must not repeat the chromosome";
			bases[i]=reverse[0]?complement(genome[len-1-at]):genome[at];
		}
		return new Window(bases, left, left+span-1, start, end);
	}

	/** Add real source context around a frozen window without moving its anchors. */
	public Window extendWindow(byte[] genome,Window original,int padding){
		require(genome!=null && original!=null && padding>=0 && !mixedStrand,
			"Only an exact reconstructed sense window may acquire additional genomic context");
		final int origin=original.sourceStart-original.start,len=genome.length;
		final int left=crossesOrigin?padding:Math.min(padding,origin);
		final int right=crossesOrigin?padding:Math.min(padding,len-origin-original.bases.length);
		require(left>=0 && right>=0,"An ordinary source window must fit its physical contig before expansion");
		final byte[] bases=new byte[original.bases.length+left+right];
		int at=(origin-left)%len;if(at<0){at+=len;}
		for(int i=0;i<bases.length;i++){
			bases[i]=reverse[0]?complement(genome[len-1-at]):genome[at];
			if(++at==len){at=0;}
		}
		for(int i=0;i<original.bases.length;i++){require(bases[left+i]==original.bases[i],"Context expansion must preserve every original window base and its orientation");}
		// Explicit origin-crossing geometry permits repeated circular context only
		// in these extra flanks. The original locus and raw anchors occur once.
		return new Window(bases,original.start+left,original.stop+left,original.sourceStart,original.sourceStop);
	}

	/** Caller's bestOffset adds this value to a sense-oriented raw anchor. */
	public static int correction(int truth, int raw){
		require(truth>=0 && raw>=0, "Endpoint correction requires valid sense coordinates");return truth-raw;
	}
	private static byte complement(byte b){
		require(b>=0 && b<AminoAcid.baseToComplementExtended.length, "Non-ASCII genomic base cannot be complemented");return AminoAcid.baseToComplementExtended[b];
	}
	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}
	public static final class Window{
		Window(byte[] b, int s, int e, int gs, int ge){bases=b;start=s;stop=e;sourceStart=gs;sourceStop=ge;}
		public final byte[] bases;
		public final int start,stop,sourceStart,sourceStop;
	}
	public final boolean mixedStrand,crossesOrigin;
	public final int[] lo,hi;
	public final boolean[] reverse;
}
