package test.assemble;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

import assemble.Contig;
import assemble.HollowBridgeJoiner;
import assemble.HollowTipConnections;
import assemble.HollowTipConnections.Connection;

/** Tiny in-memory evidence fixtures for endpoint-indexed connection discovery.
 * @author Brian Bushnell, Fischl */
public final class HollowTipConnectionsTest {

	public static void main(final String[] args){
		boolean assertions=false;
		assert(assertions=true);
		if(!assertions){throw new IllegalStateException("Connection fixtures require -ea.");}
		testOrientationAndGap();
		testMaskedInteriorAndReadBoundaries();
		testUnmappedCompetitor();
		testEndConflicts();
		testInternalCopy();
		testSliding();
		testSlidingMissingEvidence();
		testSlidingGeometry();
		System.err.println("HOLLOW_TIP_CONNECTIONS_PASS checks="+checks);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Fixtures           ----------------*/
	/*--------------------------------------------------------------*/

	/** Checks all strand combinations, positive/zero/negative gaps and noninflated support. */
	private static void testOrientationAndGap(){
		final byte[] reference=sequence(750, 739);
		for(int gap : new int[]{40, 0, -20}){
			for(int orientation=0; orientation<4; orientation++){
				final boolean leftReverse=(orientation&1)!=0, rightReverse=(orientation&2)!=0;
				byte[] left=Arrays.copyOfRange(reference, 0, 300);
				byte[] right=Arrays.copyOfRange(reference, 300+gap, 600+gap);
				if(leftReverse){left=reverse(left);}
				if(rightReverse){right=reverse(right);}
				final HollowTipConnections finder=new HollowTipConnections(contigs(left, right), 128, 128, 3);
				finder.addRead(reference);
				finder.addRead(reverse(reference));
				finder.addRead(reference);
				final ArrayList<Connection> links=finder.discover();
				check(links.size()==1, "Unique reference-supported endpoints must yield one physical connection.");
				final Connection link=links.get(0);
				check(link.leftEnd==(leftReverse ? 0 : 1) && link.rightEnd==(rightReverse ? 3 : 2),
						"Physical endpoint identities must follow the observed contig strands.");
				check(link.gap==gap && !link.conflicted(), "Signed gap must equal the fixture separation without conflicts.");
				check(link.support()==3 && link.anchorPairs()>1, "Correlated anchor windows must not inflate pair support.");
				check(finder.reads()==3 && finder.readBases()==2250 && finder.windows()==3L*(750-128+1),
						"Read/window counters must count each observed span once, regardless of strand.");
				check(finder.distinctPairs()>0 && finder.retainedWindows()<=finder.windows(), "Tip filtering must retain a valid evidence subset.");
				boolean rejected=false;
				try{finder.addRead(reference);}catch(IllegalStateException e){rejected=true;}
				check(rejected, "Discovery must freeze its evidence.");
			}
		}
	}

	/** Masked substitutions preserve endpoints; undefined spans and cross-read windows do not count. */
	private static void testMaskedInteriorAndReadBoundaries(){
		final byte[] reference=sequence(700, 751);
		final ArrayList<Contig> assembly=contigs(Arrays.copyOfRange(reference, 0, 300), Arrays.copyOfRange(reference, 340, 640));
		final byte[] spanning=Arrays.copyOfRange(reference, 244, 372);
		final HollowTipConnections masked=new HollowTipConnections(assembly, 128, 128, 3);
		for(int i=0; i<3; i++){
			final byte[] read=spanning.clone();
			read[65+i]=(byte)(read[65+i]=='A' ? 'C' : 'A');
			masked.addRead(read);
		}
		final ArrayList<Connection> links=masked.discover();
		check(links.size()==1 && links.get(0).gap==40 && links.get(0).support()==3,
				"Different masked interiors must share their literal outer-anchor support.");
		check(links.get(0).anchorPairs()==1, "An exact-span read contributes one pair, not a fabricated traversal.");
		final HollowTipConnections undefined=new HollowTipConnections(assembly, 128, 128, 1);
		final byte[] unknown=spanning.clone();
		unknown[64]='N';
		undefined.addRead(unknown);
		check(undefined.discover().isEmpty() && undefined.windows()==0,
				"This prototype deliberately excludes spans with undefined interior bases.");
		final HollowTipConnections separated=new HollowTipConnections(assembly, 128, 128, 1);
		separated.addRead(Arrays.copyOfRange(spanning, 0, 64));
		separated.addRead(Arrays.copyOfRange(spanning, 64, 128));
		check(separated.discover().isEmpty() && separated.windows()==0, "No evidence window may cross a read boundary.");
	}

	/** A supported partner absent from all tips still competes with a seemingly unique join. */
	private static void testUnmappedCompetitor(){
		final byte[] reference=sequence(700, 761);
		final ArrayList<Contig> assembly=contigs(Arrays.copyOfRange(reference, 0, 300), Arrays.copyOfRange(reference, 340, 640));
		final byte[] supported=Arrays.copyOfRange(reference, 244, 372);
		final byte[] alternate=supported.clone(), other=sequence(32, 769);
		System.arraycopy(other, 0, alternate, 96, 32);
		final HollowTipConnections finder=new HollowTipConnections(assembly, 128, 128, 3);
		for(int i=0; i<3; i++){finder.addRead(supported); finder.addRead(alternate);}
		check(finder.distinctPairs()==2, "Unmapped partners of a unique source anchor must remain in the evidence table.");
		check(finder.discover().isEmpty(), "Supported unmapped competitor must veto the apparent unique tip partner.");
		final HollowTipConnections weak=new HollowTipConnections(assembly, 128, 128, 3);
		weak.addRead(supported);
		weak.addRead(supported);
		check(weak.discover().isEmpty(), "Pair support below the explicit threshold is not a candidate.");
	}

	/** Different source anchors can suggest different neighbors; end-level conflicts must survive aggregation. */
	private static void testEndConflicts(){
		final byte[] left=sequence(400, 773), right=sequence(300, 787), other=sequence(300, 797);
		final byte[] middle=sequence(40, 809);
		final byte[] read1=join(Arrays.copyOfRange(left, 248, 400), middle, Arrays.copyOfRange(right, 0, 32));
		final byte[] read2=join(Arrays.copyOfRange(left, 328, 400), middle, Arrays.copyOfRange(other, 0, 112));
		check(read1.length==224 && read2.length==224, "Competing fixtures must contribute one disjoint anchor pair each.");
		final HollowTipConnections finder=new HollowTipConnections(contigs(left, right, other), 224, 224, 3);
		for(int i=0; i<3; i++){finder.addRead(read1); finder.addRead(read2);}
		final ArrayList<Connection> links=finder.discover();
		check(links.size()==2, "Distinct individually unambiguous pairs must remain distinct end candidates.");
		for(Connection link : links){
			check(link.leftEnd==1 && link.gap==40 && link.support()==3 && link.conflicted(),
					"An end with two supported neighbors must not pick one by iteration order.");
		}
	}

	/** Duplicated destination sequence anywhere in the assembly prevents unique endpoint assignment. */
	private static void testInternalCopy(){
		final byte[] reference=sequence(700, 811);
		final byte[] left=Arrays.copyOfRange(reference, 0, 300), right=Arrays.copyOfRange(reference, 340, 640);
		final byte[] decoy=join(sequence(100, 821), right, sequence(100, 823));
		final HollowTipConnections finder=new HollowTipConnections(contigs(left, right, decoy), 128, 128, 3);
		final byte[] spanning=Arrays.copyOfRange(reference, 244, 372);
		for(int i=0; i<3; i++){finder.addRead(spanning);}
		check(finder.discover().isEmpty(), "An internal reverse/forward copy is not a unique destination tip.");
	}

	/** Inland seeding alone misses offset-only branches; sliding must reject them on either strand. */
	private static void testSliding(){
		final byte[] reference=sequence(750, 829);
		for(int orientation=0; orientation<4; orientation++){
			byte[] left=Arrays.copyOfRange(reference, 0, 300), right=Arrays.copyOfRange(reference, 340, 640);
			if((orientation&1)!=0){left=reverse(left);}
			if((orientation&2)!=0){right=reverse(right);}
			for(int competitor=0; competitor<3; competitor++){
				final ArrayList<Contig> assembly=contigs(left, right);
				final HollowTipConnections legacy=new HollowTipConnections(assembly, 128, 128, 3, 12, 12, 0);
				final HollowTipConnections sliding=new HollowTipConnections(assembly, 128, 128, 3, 12, 12, 10);
				final byte[] alternate=Arrays.copyOfRange(reference, competitor==2 ? 251 : 261, competitor==2 ? 379 : 389);
				final int changed=competitor==2 ? 0 : 127;
				alternate[changed]=(byte)(alternate[changed]=='A' ? 'C' : 'A');
				for(int i=0; i<3; i++){
					legacy.addRead(i==1 ? reverse(reference) : reference);
					sliding.addRead(i==1 ? reverse(reference) : reference);
					if(competitor>0){legacy.addRead(alternate); sliding.addRead(alternate);}
				}
				final ArrayList<Connection> old=legacy.discover(), links=sliding.discover();
				check(old.size()==1 && links.size()==1, "Only the exact12bp-inland seed may define the fixture connection.");
				final Connection c=links.get(0);
				check(c.anchorPairs()==1 && c.gap==40 && !c.conflicted(), "Inland selection must preserve signed geometry.");
				check(c.slideChecked()==21 && c.slideMinCount()==3, "Every offset, including both endpoints, needs expected support.");
				check(competitor==0 ? c.slideReason()==null : "sliding_branch".equals(c.slideReason()),
						"Offset-only incoming/outgoing alternatives must veto without hiding the seed candidate.");
				final HollowBridgeJoiner joiner=new HollowBridgeJoiner(assembly, sliding, links, 3, 4);
				for(int i=0; i<3; i++){joiner.addRead(reference);}
				joiner.join();
				check(joiner.joins()==(competitor==0 ? 1 : 0), "Sliding veto must reach the actual assembly consumer.");
			}
		}
		// Both shifted arms are duplicated internally; they cannot seed, but their partners remain negative evidence.
		final byte[] decoy=sequence(600, 839);
		System.arraycopy(reference, 261, decoy, 200, 32);
		System.arraycopy(reference, 357, decoy, 400, 32);
		final ArrayList<Contig> assembly=contigs(Arrays.copyOfRange(reference, 0, 300), Arrays.copyOfRange(reference, 340, 640), decoy);
		final HollowTipConnections finder=new HollowTipConnections(assembly, 128, 128, 3, 12, 12, 10);
		final byte[] alternate=Arrays.copyOfRange(reference, 261, 389);
		alternate[127]=(byte)(alternate[127]=='A' ? 'C' : 'A');
		for(int i=0; i<3; i++){finder.addRead(reference); finder.addRead(alternate);}
		final ArrayList<Connection> links=finder.discover();
		check(links.size()==1 && "sliding_branch".equals(links.get(0).slideReason()),
				"Nonunique shifted arms must not disappear from the branch evidence index.");
	}

	/** A center-only observation cannot establish sliding support; boundary seeds cannot pretend to slide. */
	private static void testSlidingMissingEvidence(){
		final byte[] reference=sequence(750, 853);
		final ArrayList<Contig> assembly=contigs(Arrays.copyOfRange(reference, 0, 300), Arrays.copyOfRange(reference, 340, 640));
		final HollowTipConnections finder=new HollowTipConnections(assembly, 128, 128, 3, 12, 12, 10);
		for(int i=0; i<3; i++){finder.addRead(Arrays.copyOfRange(reference, 256, 384));}
		final ArrayList<Connection> links=finder.discover();
		check(links.size()==1 && "sliding_low_support".equals(links.get(0).slideReason()),
				"Missing neighboring spans must reject with low support, not count as unbranched proof.");
		check(links.get(0).slideChecked()==21 && links.get(0).slideMinCount()==0,
				"Missing evidence must remain visible in the offset diagnostics.");
		final HollowTipConnections edge=new HollowTipConnections(assembly, 104, 128, 3, 0, 0, 10);
		for(int i=0; i<3; i++){edge.addRead(reference);}
		final ArrayList<Connection> edges=edge.discover();
		check(edges.size()==1 && "sliding_margin".equals(edges.get(0).slideReason()),
				"Both arms at their tips cannot supply a symmetric known-sequence sliding window.");
		boolean rejected=false;
		try{new HollowTipConnections(assembly, 128, 128, 3, 20, 10, 10);}
		catch(IllegalArgumentException e){rejected=true;}
		check(rejected, "An inverted inland interval is not a usable seed band.");
	}

	/** Sliding preserves positive gaps, abutments and overlaps at all four endpoint orientations. */
	private static void testSlidingGeometry(){
		final byte[] reference=sequence(750, 857);
		for(int gap : new int[]{40, 0, -20}){
			for(int orientation=0; orientation<4; orientation++){
				byte[] left=Arrays.copyOfRange(reference, 0, 300), right=Arrays.copyOfRange(reference, 300+gap, 600+gap);
				if((orientation&1)!=0){left=reverse(left);}
				if((orientation&2)!=0){right=reverse(right);}
				final HollowTipConnections finder=new HollowTipConnections(contigs(left, right), 88+gap, 128, 3, 12, 12, 10);
				for(int i=0; i<3; i++){finder.addRead(reference);}
				final ArrayList<Connection> links=finder.discover();
				check(links.size()==1 && links.get(0).gap==gap && links.get(0).slideReason()==null,
						"Signed overlap/gap and strand must survive sliding verification.");
				check(links.get(0).slideChecked()==21, "Every valid geometry must check both full10bp flanks.");
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Test helpers         ----------------*/
	/*--------------------------------------------------------------*/

	private static ArrayList<Contig> contigs(final byte[]... sequences){
		final ArrayList<Contig> contigs=new ArrayList<Contig>();
		for(byte[] sequence : sequences){contigs.add(new Contig(sequence, contigs.size()));}
		return contigs;
	}

	/** Small deterministic in-memory fixtures, not new benchmark libraries. */
	private static byte[] sequence(final int length, final int seed){
		final byte[] bases=new byte[length], alphabet={'A','C','G','T'};
		final Random random=new Random(seed);
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

	private static byte[] join(final byte[] a, final byte[] b, final byte[] c){
		final byte[] result=new byte[a.length+b.length+c.length];
		System.arraycopy(a, 0, result, 0, a.length);
		System.arraycopy(b, 0, result, a.length, b.length);
		System.arraycopy(c, 0, result, a.length+b.length, c.length);
		return result;
	}

	private static void check(final boolean condition, final String message){
		checks++;
		if(!condition){throw new AssertionError(message);}
	}

	private static int checks;
}
