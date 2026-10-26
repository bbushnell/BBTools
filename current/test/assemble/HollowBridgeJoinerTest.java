package test.assemble;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

import assemble.Contig;
import assemble.HollowBridgeJoiner;
import assemble.HollowTipConnections;
import assemble.HollowTipConnections.Connection;

/** Tiny in-memory tests for read-backed gap reconstruction and oriented path spelling.
 * @author Brian Bushnell, Fischl */
public final class HollowBridgeJoinerTest {

	public static void main(final String[] args){
		boolean assertions=false;
		assert(assertions=true);
		if(!assertions){throw new IllegalStateException("Join fixtures require -ea.");}
		testOrientations();
		testScaffolds();
		testAmbiguityAndFlanks();
		testRepeatedPair();
		testChain();
		testCycle();
		testOuterFlanks();
		testFlankDistance();
		System.err.println("HOLLOW_BRIDGE_JOINER_PASS checks="+checks);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Fixtures           ----------------*/
	/*--------------------------------------------------------------*/

	/** Unknown-gap spelling retains both complete contigs, even when the inferred gap is negative. */
	private static void testScaffolds(){
		final byte[] reference=sequence(1000, 971);
		for(int gap : new int[]{40, 0, -20}){
			for(int orientation=0; orientation<4; orientation++){
				final byte[] left=Arrays.copyOfRange(reference, 0, 300);
				final byte[] right=Arrays.copyOfRange(reference, 300+gap, 600+gap);
				final ArrayList<Contig> input=contigs((orientation&1)==0 ? left : reverse(left),
						(orientation&2)==0 ? right : reverse(right));
				final HollowTipConnections finder=find(input, reference);
				final HollowBridgeJoiner joiner=new HollowBridgeJoiner(input, finder, finder.discover(), 3, 4);
				for(int i=0; i<3; i++){joiner.addRead(i==1 ? reverse(reference) : reference);}
				final ArrayList<Contig> output=joiner.join();
				check(output.size()==1 && joiner.joins()==1, "Default spelling must retain the accepted link.");
				final byte[] expected=new byte[610];
				System.arraycopy(left, 0, expected, 0, 300);
				Arrays.fill(expected, 300, 310, (byte)'N');
				System.arraycopy(right, 0, expected, 310, 300);
				check(Arrays.equals(output.get(0).bases, expected),
						"Scaffolding must retain every input base with exactly ten unknown bases, never fill or trim.");
			}
		}
		final ArrayList<Contig> input=contigs(Arrays.copyOfRange(reference, 0, 300),
				reverse(Arrays.copyOfRange(reference, 340, 640)), Arrays.copyOfRange(reference, 680, 980));
		final HollowTipConnections finder=find(input, reference);
		final HollowBridgeJoiner joiner=new HollowBridgeJoiner(input, finder, finder.discover(), 3, 4, 128, 16);
		for(int i=0; i<3; i++){joiner.addRead(reference);}
		final ArrayList<Contig> output=joiner.join();
		check(output.size()==1 && joiner.joins()==2, "The default guarded constructor must scaffold both links.");
		final byte[] expected=new byte[920];
		for(int i=0; i<3; i++){System.arraycopy(reference, i*340, expected, i*310, 300);}
		Arrays.fill(expected, 300, 310, (byte)'N');
		Arrays.fill(expected, 610, 620, (byte)'N');
		check(Arrays.equals(output.get(0).bases, expected), "An oriented chain must preserve all three full contigs.");
	}

	/** Every strand combination must produce the same complete sequence across gaps and overlaps. */
	private static void testOrientations(){
		final byte[] reference=sequence(750, 971);
		for(int gap : new int[]{40, 0, -20}){
			for(int orientation=0; orientation<4; orientation++){
				byte[] left=Arrays.copyOfRange(reference, 0, 300), right=Arrays.copyOfRange(reference, 300+gap, 600+gap);
				if((orientation&1)!=0){left=reverse(left);}
				if((orientation&2)!=0){right=reverse(right);}
				final ArrayList<Contig> contigs=contigs(left, right);
				final HollowTipConnections finder=find(contigs, reference);
				final ArrayList<Connection> links=finder.discover();
				check(links.size()==1, "One unique fixture connection is required before reconstruction.");
				final HollowBridgeJoiner joiner=new HollowBridgeJoiner(contigs, finder, links, 3, 4, 0, 16, true);
				joiner.addRead(reference); joiner.addRead(reverse(reference)); joiner.addRead(reference);
				final ArrayList<Contig> output=joiner.join();
				check(output.size()==1 && joiner.joins()==1, "Supported reciprocal geometry must spell one path.");
				check(Arrays.equals(output.get(0).bases, Arrays.copyOfRange(reference, 0, 600+gap)),
						"Gap/overlap spelling must reconstruct the exact original reference in all contig orientations.");
				check(joiner.evidence().get(0).support()==3, "Many anchor windows per read must not multiply sequence votes.");
				final HollowBridgeJoiner guarded=new HollowBridgeJoiner(contigs, finder, links, 3, 4, 128, 16, true);
				guarded.addRead(reference); guarded.addRead(reverse(reference)); guarded.addRead(reference);
				check(guarded.join().size()==1 && guarded.evidence().get(0).outerConfirmed()==3,
						"Outer-flank checks must preserve clean joins in every contig/read orientation and gap geometry.");
			}
		}
	}

	/** Tied interiors and mismatching known flanks must not be converted into an invented consensus. */
	private static void testAmbiguityAndFlanks(){
		final byte[] reference=sequence(700, 977), read=Arrays.copyOfRange(reference, 244, 372);
		final byte[] alternate=read.clone(), badFlank=read.clone();
		mutate(alternate, 70); mutate(badFlank, 40);
		final ArrayList<Contig> contigs=contigs(Arrays.copyOfRange(reference, 0, 300), Arrays.copyOfRange(reference, 340, 640));
		final HollowTipConnections finder=find(contigs, read);
		final ArrayList<Connection> links=finder.discover();
		check(links.size()==1, "The selected exact-span pair must identify the fixture gap.");
		final HollowBridgeJoiner tied=new HollowBridgeJoiner(contigs, finder, links, 3, 4, 0, 16, true);
		for(int i=0; i<3; i++){tied.addRead(read); tied.addRead(alternate);}
		check(tied.join().size()==2 && tied.evidence().get(0).reason().equals("gap_sequence_ambiguity"),
				"Equally supported whole-gap sequences must abstain.");
		final HollowBridgeJoiner bad=new HollowBridgeJoiner(contigs, finder, links, 3, 4, 0, 16, true);
		for(int i=0; i<3; i++){bad.addRead(badFlank);}
		check(bad.join().size()==2 && bad.evidence().get(0).flankRejected()==3,
				"Exact endpoint arms do not excuse disagreement with known internal contig sequence.");
		final HollowBridgeJoiner dominant=new HollowBridgeJoiner(contigs, finder, links, 3, 4, 0, 16, true);
		for(int i=0; i<12; i++){dominant.addRead(read);}
		for(int i=0; i<3; i++){dominant.addRead(alternate);}
		check(dominant.join().size()==1 && dominant.evidence().get(0).support()==12,
				"Exact vote-ratio equality permits the supported complete gap, not a mosaic.");
	}

	/** A repeated selector in one record cannot masquerade as multiple independent observations. */
	private static void testRepeatedPair(){
		final byte[] reference=sequence(700, 983), read=Arrays.copyOfRange(reference, 244, 372);
		final byte[] twice=new byte[256];
		System.arraycopy(read, 0, twice, 0, 128); System.arraycopy(read, 0, twice, 128, 128);
		final ArrayList<Contig> contigs=contigs(Arrays.copyOfRange(reference, 0, 300), Arrays.copyOfRange(reference, 340, 640));
		final HollowTipConnections finder=find(contigs, read);
		final HollowBridgeJoiner joiner=new HollowBridgeJoiner(contigs, finder, finder.discover(), 3, 4, 0, 16, true);
		for(int i=0; i<3; i++){joiner.addRead(twice);}
		check(joiner.join().size()==2 && joiner.evidence().get(0).support()==0,
				"A second selected-pair occurrence revokes the first occurrence's vote from that read record.");
	}

	/** A short middle contig may participate in two valid joins within one linear path. */
	private static void testChain(){
		final byte[] reference=sequence(1000, 991);
		final ArrayList<Contig> contigs=contigs(Arrays.copyOfRange(reference, 0, 300),
				reverse(Arrays.copyOfRange(reference, 340, 640)), Arrays.copyOfRange(reference, 680, 980));
		final HollowTipConnections finder=find(contigs, reference);
		final ArrayList<Connection> links=finder.discover();
		check(links.size()==2, "Both sides of the middle contig must be discovered.");
		final HollowBridgeJoiner joiner=new HollowBridgeJoiner(contigs, finder, links, 3, 4, 0, 16, true);
		for(int i=0; i<3; i++){joiner.addRead(reference);}
		final ArrayList<Contig> output=joiner.join();
		check(output.size()==1 && joiner.joins()==2 && Arrays.equals(output.get(0).bases, Arrays.copyOf(reference, 980)),
				"An oriented chain must include each contig exactly once and spell both gaps correctly.");
	}

	/** A degree-two closed component is preserved rather than arbitrarily cut or circularized. */
	private static void testCycle(){
		final byte[] reference=sequence(960, 997);
		final ArrayList<Contig> contigs=contigs(Arrays.copyOfRange(reference, 0, 300),
				Arrays.copyOfRange(reference, 320, 620), Arrays.copyOfRange(reference, 640, 940));
		final byte[][] reads=new byte[3][128];
		for(int r=0; r<3; r++){
			for(int i=0; i<128; i++){reads[r][i]=reference[(240+320*r+i)%960];}
		}
		final HollowTipConnections finder=new HollowTipConnections(contigs, 128, 128, 3);
		for(byte[] read : reads){for(int i=0; i<3; i++){finder.addRead(read);}}
		final ArrayList<Connection> links=finder.discover();
		check(links.size()==3, "Circular fixture must expose three different physical connections.");
		final HollowBridgeJoiner joiner=new HollowBridgeJoiner(contigs, finder, links, 3, 4, 0, 16, true);
		for(byte[] read : reads){for(int i=0; i<3; i++){joiner.addRead(read);}}
		final ArrayList<Contig> output=joiner.join();
		check(output.size()==3 && joiner.joins()==0, "Cycles are not supported by this conservative linear-path prototype.");
		for(int i=0; i<3; i++){check(output.get(i)==contigs.get(i), "An unjoined contig must retain its original object and sequence.");}
	}

	/** Conflicts are candidate-level evidence; single edits, Ns, clipping and repeated selectors stay conservative. */
	private static void testOuterFlanks(){
		final byte[] reference=sequence(700, 1009), seed=Arrays.copyOfRange(reference, 244, 372);
		final byte[] read=Arrays.copyOfRange(reference, 200, 420), bad=read.clone(), unknown=read.clone();
		mutate(bad, 20); mutate(bad, 25);
		Arrays.fill(unknown, 0, 44, (byte)'N');
		Arrays.fill(unknown, 172, unknown.length, (byte)'N');
		final ArrayList<Contig> contigs=contigs(Arrays.copyOfRange(reference, 0, 300), Arrays.copyOfRange(reference, 340, 640));
		final HollowTipConnections finder=find(contigs, seed);
		final ArrayList<Connection> links=finder.discover();
		check(links.size()==1, "A fixed short seed isolates outside-window evidence from discovery.");
		for(boolean reverse : new boolean[]{false, true}){
			final HollowBridgeJoiner guarded=new HollowBridgeJoiner(contigs, finder, links, 3, 4, 128, 16, true);
			for(int i=0; i<3; i++){
				guarded.addRead(reverse ? reverse(bad) : bad);
				guarded.addRead(reverse ? reverse(read) : read);
			}
			check(guarded.join().size()==2 && guarded.evidence().get(0).reason().equals("outer_flank_conflict"),
					"Three contradictory reads veto the candidate even with three valid gap votes.");
			check(guarded.evidence().get(0).outerContradictory()==3 && guarded.evidence().get(0).outerConfirmed()==3,
					"Both read orientations must retain exactly one signed vote per record.");
		}
		final HollowBridgeJoiner tolerant=new HollowBridgeJoiner(contigs, finder, links, 3, 4, 128, 16, true);
		final byte[] sub=read.clone(); mutate(sub, 20);
		final byte[] insertion=new byte[read.length+1], deletion=new byte[read.length-1];
		System.arraycopy(read, 0, insertion, 0, 20); insertion[20]='A';
		System.arraycopy(read, 20, insertion, 21, read.length-20);
		System.arraycopy(read, 0, deletion, 0, 20);
		System.arraycopy(read, 21, deletion, 20, read.length-21);
		tolerant.addRead(sub); tolerant.addRead(reverse(insertion)); tolerant.addRead(deletion);
		tolerant.addRead(unknown); tolerant.addRead(seed);
		check(tolerant.join().size()==1 && tolerant.evidence().get(0).outerConfirmed()==3
				&& tolerant.evidence().get(0).outerContradictory()==0 && tolerant.evidence().get(0).support()==5,
				"One edit per flank is tolerated; Ns and clipped flanks are not falsely counted as confirmed or conflicting.");
		final HollowBridgeJoiner ratio=new HollowBridgeJoiner(contigs, finder, links, 3, 4, 128, 16, true);
		for(int i=0; i<12; i++){ratio.addRead(read);}
		for(int i=0; i<3; i++){ratio.addRead(bad);}
		check(ratio.join().size()==1 && ratio.evidence().get(0).support()==12,
				"Exact confirmed/conflicting ratio equality accepts, while conflicting records never vote for a gap.");
		final HollowBridgeJoiner duplicate=new HollowBridgeJoiner(contigs, finder, links, 3, 4, 128, 16, true);
		final byte[] twice=new byte[2*bad.length];
		System.arraycopy(bad, 0, twice, 0, bad.length); System.arraycopy(bad, 0, twice, bad.length, bad.length);
		for(int i=0; i<3; i++){duplicate.addRead(twice); duplicate.addRead(read);}
		check(duplicate.join().size()==1 && duplicate.evidence().get(0).outerContradictory()==0,
				"Repeated selected pairs revoke a negative vote as well as a positive gap vote.");
	}

	/** Exhaustive binary words compare the allocation-free one-edit test with an independent full DP oracle. */
	private static void testFlankDistance(){
		for(int length=1; length<=5; length++){
			for(int q=0; q<(1<<length); q++){
				for(int t=0; t<(1<<(length+1)); t++){
					final byte[] query=new byte[length], target=new byte[length+1];
					for(int i=0; i<length; i++){query[i]=(byte)((q>>>i)&1);}
					for(int i=0; i<=length; i++){target[i]=(byte)((t>>>i)&1);}
					final int[][] dp=new int[length+1][length+2];
					for(int i=0; i<=length; i++){dp[i][0]=i;}
					for(int j=0; j<=length+1; j++){dp[0][j]=j;}
					for(int i=1; i<=length; i++){
						for(int j=1; j<=length+1; j++){
							dp[i][j]=Math.min(dp[i-1][j-1]+(query[i-1]==target[j-1] ? 0 : 1),
									Math.min(dp[i-1][j]+1, dp[i][j-1]+1));
						}
					}
					int best=length;
					for(int j=0; j<=length+1; j++){best=Math.min(best, dp[length][j]);}
					check(HollowBridgeJoiner.flankAgrees(query, length, target, length+1)==(best<=1),
							"Single-edit shortcut must equal anchored full-query/free-reference-end edit distance.");
				}
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Test helpers         ----------------*/
	/*--------------------------------------------------------------*/

	private static HollowTipConnections find(final ArrayList<Contig> contigs, final byte[] read){
		final HollowTipConnections finder=new HollowTipConnections(contigs, 128, 128, 3);
		for(int i=0; i<3; i++){finder.addRead(read);}
		return finder;
	}

	private static ArrayList<Contig> contigs(final byte[]... sequences){
		final ArrayList<Contig> result=new ArrayList<Contig>();
		for(byte[] sequence : sequences){result.add(new Contig(sequence, result.size()));}
		return result;
	}

	/** Tiny deterministic in-memory fixtures, not additional benchmark libraries. */
	private static byte[] sequence(final int length, final int seed){
		final Random random=new Random(seed);
		final byte[] bases=new byte[length], alphabet={'A','C','G','T'};
		for(int i=0; i<length; i++){bases[i]=alphabet[random.nextInt(4)];}
		return bases;
	}

	private static byte[] reverse(final byte[] bases){
		final byte[] result=new byte[bases.length];
		for(int i=0; i<bases.length; i++){
			final byte b=bases[bases.length-1-i];
			result[i]=(byte)(b=='A' ? 'T' : b=='C' ? 'G' : b=='G' ? 'C' : 'A');
		}
		return result;
	}

	private static void mutate(final byte[] bases, final int position){bases[position]=(byte)(bases[position]=='A' ? 'C' : 'A');}
	private static void check(final boolean condition, final String message){
		checks++;
		if(!condition){throw new AssertionError(message);}
	}
	private static int checks;
}
