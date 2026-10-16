package assemble;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.ObjectIntMap;
import parse.LineParser1;
import structures.ByteBuilder;

/** Independent small-fixture count reconstruction for live branch traces. @author Fischl */
public final class LiveBranchObserverTest {

	/** Reads a finished native trace after its read count table has been released. */
	public static void main(final String[] args){
		if(args.length==3 && args[0].equals("generate")){generate(args[1], args[2]); return;}
		if(args.length!=3){throw new IllegalArgumentException("Expected reads.fa trace.tsv K; use controlevery=1 branchmult1=1000 branchlower=0 mce=2 wash=f.");}
		checkTrace(args[0], args[1], Integer.parseInt(args[2]));
	}

	/** Uses 300bp reads so a central substitution has complete K95 support on both sides. */
	private static void generate(final String reference, final String output){
		final String genome=FusionJoinDiagnostic.sequences(reference).get(0);
		if(genome.length()<300){throw new IllegalArgumentException("Branch fixture reference must contain a300bp read.");}
		final ByteStreamWriter writer=new ByteStreamWriter(output, false, false, false);
		final ByteBuilder row=new ByteBuilder();
		writer.start();
		int id=0;
		try{
			for(int pos=0; pos+300<=genome.length(); pos+=10){
				for(int copy=0; copy<2; copy++){
					row.clear().append('>').append(id++).nl();
					for(int i=0; i<300; i++){
						char base=genome.charAt(pos+i);
						if(copy==0 && pos%70==0 && i==149){base=base=='A' ? 'C' : 'A';}
						row.append(base);
					}
					writer.print(row.nl());
				}
			}
		}finally{if(writer.poisonAndWait()){throw new IllegalStateException("Failed to write300bp branch fixture.");}}
		System.err.println("LIVE_BRANCH_READ_FIXTURE reads="+id+" length=300 central_substitution=true");
	}

	/** Reconstructs every observed count from literal input substrings, not Tadpole APIs. */
	private static void checkTrace(final String reads, final String trace, final int k){
		if(k<24 || k>127){throw new IllegalArgumentException("Fixture K must lie in24–127.");}
		final ObjectIntMap<String> expected=new ObjectIntMap<String>(String.class);
		final ArrayList<String> sequences=FusionJoinDiagnostic.sequences(reads);
		for(int r=0; r<sequences.size(); r+=2){
			final String sequence=sequences.get(r);
			for(int p=0; p+k<=sequence.length(); p++){
				final String word=canonical(sequence.substring(p, p+k));
				// Explicit absent-key handling: ObjectIntMap.increment currently adds
				// to resize()'s -1 sentinel on first insertion, rather than to zero.
				expected.put(word, Math.max(0, expected.get(word))+1);
			}
		}
		long observations=0, squares=0;
		for(int count : expected.values()){
			if(count>0){observations+=count; squares+=(long)count*count;}
		}
		final ByteFile input=ByteFile.makeByteFile(trace, true);
		final LineParser1 fields=new LineParser1('\t');
		long rows=0, branchRows=0, interior=0;
		boolean complete=false;
		try{
			check("#livebranch-observation-v2".equals(text(input.nextLine())), "Missing trace schema.");
			check(LiveBranchObserver.HEADER.equals(text(input.nextLine())), "Wrong trace header.");
			for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
				check(!complete, "Data appeared after terminal evidence.");
				fields.set(line);
				if(fields.termEquals("#complete", 0)){
					check(fields.terms()==9 && fields.parseLong(1)==rows, "Completion must cover every trace row.");
					check(fields.parseLong(2)==expected.size(), "Live table unique count differs from literal input.");
					// The two terminal doubles use round-trip text, which may contain E.
					// Native fast parseDouble intentionally excludes exponent notation.
					check(Double.parseDouble(fields.parseString(3))==observations &&
							Double.parseDouble(fields.parseString(4))==squares, "Global depth moments differ from literal input.");
					check(fields.parseInt(5)==1, "Fixture must collect all native controls.");
					check(fields.parseInt(6)==k && fields.parseInt(7)==2 && fields.parseDouble(8)==0.5,
							"Trace must bind its K, mce and hard-stop screen explicitly.");
					complete=true;
					continue;
				}
				check(fields.terms()==19 && fields.parseInt(1)==k, "Wrong trace width/K.");
				final String word=fields.parseString(0);
				check(word.length()==k && word.equals(canonical(word)), "Trace word is not canonical.");
				check(fields.parseInt(2)==count(expected, word), "Wrong node count for "+word+
						": observed="+fields.parseInt(2)+", literal="+count(expected, word));
				final int[] left=new int[4], right=new int[4];
				for(int base=0; base<4; base++){
					left[base]=count(expected, "ACGT".charAt(base)+word.substring(0, k-1));
					right[base]=count(expected, word.substring(1)+"ACGT".charAt(base));
					check(fields.parseInt(3+base)==left[base] && fields.parseInt(7+base)==right[base], "Strand/base-indexed neighbor count differs.");
				}
				final boolean lb=branch(left), rb=branch(right);
				check(fields.parseInt(11)==(lb ? 1 : 0) && fields.parseInt(12)==(rb ? 1 : 0), "Native branch flags differ after canonicalization.");
				final int arrivalSide=fields.parseInt(13), arrival=fields.parseInt(14), entry=fields.parseInt(15);
				check(entry==0 || entry==1, "Entry provenance must be binary.");
				if(entry==1){check(arrivalSide==-1 && arrival==-1, "Entry must not invent a predecessor.");}
				else{
					check((arrivalSide==0 || arrivalSide==1) && arrival>=0 && arrival<4, "Interior traversal lost its predecessor.");
					check((arrivalSide==0 ? left[arrival] : right[arrival])>=2, "Recorded incoming path has insufficient read support.");
					interior++;
				}
				check(fields.parseInt(17)==LiveBranchPolicy.level(fields.parseInt(2), left, right, lb, rb, 2, 0.5), "Recorded screen differs from symmetric policy.");
				check(fields.parseInt(18)==1, "All-control fixture must have sampling denominator1.");
				rows++;
				if(lb || rb){branchRows++;}
			}
		}finally{if(input.close()){throw new IllegalStateException("Failed to read native observer trace.");}}
		check(complete && rows>0 && branchRows>0 && interior>0, "Fixture must exercise completed branch/interior observations: complete="+
				complete+", rows="+rows+", branch="+branchRows+", interior="+interior);
		System.err.println("LIVE_BRANCH_OBSERVER_PASS k="+k+" rows="+rows+" branch_rows="+branchRows+" interior="+interior);
	}

	/** Matches only the explicitly fixed native fixture thresholds, independently of Tadpole. */
	private static boolean branch(final int[] counts){
		assert(counts.length==4) : "Fixture branch truth uses four nucleotide depths.";
		int max=0, second=0;
		for(int count : counts){
			if(count>max){second=max; max=count;}
			else if(count>second){second=count;}
		}
		return second>0 && second*1000L>=max;
	}

	/** Absent literal words have zero support, not the map's -1 sentinel. */
	private static int count(final ObjectIntMap<String> map, final String word){
		assert(!word.isEmpty()) : "Only complete nonempty Kmers belong in count fixtures.";
		return Math.max(0, map.get(canonical(word)));
	}

	/** Lexical orientation is independent of packed-word layout and hash functions. */
	private static String canonical(final String word){
		assert(!word.isEmpty()) : "Canonicalization requires a complete observed word.";
		final StringBuilder reverse=new StringBuilder(word.length());
		for(int i=word.length()-1; i>=0; i--){reverse.append("TGCA".charAt("ACGT".indexOf(word.charAt(i))));}
		final String rc=reverse.toString();
		return word.compareTo(rc)<=0 ? word : rc;
	}

	/** Header-only conversion; raw data rows use the native line parser. */
	private static String text(final byte[] bytes){return bytes==null ? null : new String(bytes, StandardCharsets.US_ASCII);}

	/** Fails loudly without relying on a disabled assert for fixture verdicts. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
}
