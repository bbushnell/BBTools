package prok;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import dna.AminoAcid;
import map.LongIntMap;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Tools;
import stream.Read;
import structures.ByteBuilder;
import structures.IntList;
import structures.LongList;

/** Counts original occurrences and distinct gene-copy support of an ordered
 * seed extension. First-two selection credits expose redundant batch choices;
 * genes and keys are never removed, reordered or reselected by this report.
 * @author Brian Bushnell, Diona
 */
public final class ReportNcrnaSeedSupport{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, ReportNcrnaSeedSupport.class, false);
		try{
			new ReportNcrnaSeedSupport(pp.args).process();
		}finally{
			Shared.closeStream(pp.outstream);
		}
	}

	private ReportNcrnaSeedSupport(final String[] args){
		Shared.setThreads(1);
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(String arg:args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			final String a=pair[0].toLowerCase(java.util.Locale.ROOT);
			final String b=pair.length>1 ? pair[1] : null;
			if(a.equals("seedset")){initialPath=b;}
			else if(a.equals("genes")){genes=b.split(",", -1);}
			else if(a.equals("counts")){geneCounts=Parse.parseIntArray(b, ",");}
			else if(a.equals("expectedinitial")){expectedInitial=NcrnaNumbers.intValue(a, b);}
			else if(a.equals("expectedkeys")){expectedKeys=NcrnaNumbers.intValue(a, b);}
			else{require(parser.parse(arg, a, b), "Unknown support-report argument: "+arg);}
		}
		seedPath=parser.in1;
		out=parser.out1;
		require(seedPath!=null && initialPath!=null && out!=null && parser.in2==null && parser.out2==null
			&& !parser.overwrite && !parser.append && parser.maxReads<0 && parser.samplerate==1 && Shared.threads()==1,
			"Supply the complete ordered in=seeds seedset=fixed_prefix genes=originals counts=... out=fresh_prefix t=1");
		require(expectedInitial>0 && expectedKeys>=expectedInitial && genes!=null && geneCounts!=null
			&& genes.length==geneCounts.length && genes.length>0, "Declare exact seed and original-gene populations");
		for(int n:geneCounts){require(n>0, "Every original-gene file must have a positive declared count");}
		require(Tools.testInputFiles(false, true, seedPath, initialPath) && Tools.testInputFiles(false, true, genes),
			"All ordered keys and original gene copies must be readable");
		RrnaResourceIO.fresh(out, ".keys.tsv", ".summary.tsv");
		Read.TO_UPPER_CASE=true;
	}

	private void process(){
		final LongList initial=loadKeys(initialPath, expectedInitial);
		final LongList keys=loadKeys(seedPath, expectedKeys);
		for(int i=0; i<initial.size; i++){
			require(initial.get(i)==keys.get(i), "Every fixed initial key and its order must be preserved");
		}
		final LongIntMap index=new LongIntMap(expectedKeys*2);
		for(int i=0; i<keys.size; i++){index.put(keys.get(i), i);}
		final long[] occurrences=new long[expectedKeys];
		final long[] copies=new long[expectedKeys];
		final long[] credits=new long[expectedKeys];
		final long[][] coverage=new long[2][3];
		final IntList matched=new IntList();
		long total=0;
		for(int file=0; file<genes.length; file++){
			final long[] seen={0};
			final Throwable[] failure={null};
			RrnaResourceIO.stream(genes[file], read->{
				if(failure[0]!=null){return;}
				try{
					matched.clear();
					long word=0;
					int valid=0;
					for(byte base:read.bases){
						final int code=AminoAcid.baseToNumber[base&255];
						if(code<0){
							word=0;
							valid=0;
							continue;
						}
						word=((word<<2)|code)&MASK;
						if(++valid<17){continue;}
						final int selected=index.get(word);
						if(selected>=0){
							occurrences[selected]++;
							matched.add(selected);
						}
					}
					Arrays.sort(matched.array, 0, matched.size);
					int previous=-1;
					int distinct=0;
					int initialHits=0;
					for(int i=0; i<matched.size; i++){
						final int selected=matched.array[i];
						if(selected==previous){continue;}
						previous=selected;
						copies[selected]++;
						distinct++;
						if(selected<expectedInitial){initialHits++;}
						if(distinct<=2){credits[selected]++;}
					}
					coverage[0][Math.min(2, initialHits)]++;
					coverage[1][Math.min(2, distinct)]++;
					seen[0]++;
				}catch(RuntimeException | Error error){
					failure[0]=error;
				}
			});
			rethrow(failure[0]);
			require(seen[0]==geneCounts[file], "Original gene-copy denominator changed: "+genes[file]);
			total+=seen[0];
		}
		long shared=0;
		long redundant=0;
		long filled=0;
		try(Output output=new Output(out+".keys.tsv")){
			output.writer.println("selection_rank1\tkey\tfixed_initial\toriginal_occurrences\tgene_copies\tfirst_two_slot_credit");
			for(int i=0; i<keys.size; i++){
				final boolean fixed=i<expectedInitial;
				require(fixed || occurrences[i]>0, "An added key absent from the original genes cannot be an extension choice");
				require(occurrences[i]>=copies[i] && copies[i]>=credits[i], "Occurrence, copy and coverage-credit counts must be nested");
				if(!fixed){
					if(copies[i]>1){shared++;}
					if(credits[i]==0){redundant++;}
					filled+=credits[i];
				}
				output.writer.print(row.clear().append(i+1).tab().append(word(keys.get(i))).tab().append(fixed)
					.tab().append(occurrences[i]).tab().append(copies[i]).tab().append(credits[i]).nl());
			}
		}
		final long missing=2*coverage[0][0]+coverage[0][1];
		final long remaining=2*coverage[1][0]+coverage[1][1];
		require(filled==missing-remaining, "Every filled missing slot must belong to one added key in recorded selection order");
		try(Output output=new Output(out+".summary.tsv")){
			output.writer.println("genes\tinitial_keys\ttotal_keys\tadded_keys\tadded_keys_multi_copy\tadded_keys_single_copy\tadded_keys_zero_slot_credit\tinitial_missing_slots\tfilled_slots\tremaining_slots\tbefore_zero\tbefore_one\tbefore_two_plus\tafter_zero\tafter_one\tafter_two_plus\tqualification");
			row.clear().append(total).tab().append(expectedInitial).tab().append(expectedKeys)
				.tab().append(expectedKeys-expectedInitial).tab().append(shared)
				.tab().append(expectedKeys-expectedInitial-shared).tab().append(redundant)
				.tab().append(missing).tab().append(filled).tab().append(remaining);
			for(int side=0; side<2; side++){
				for(long n:coverage[side]){row.tab().append(n);}
			}
			output.writer.print(row.append("\tORIGINAL_OCCURRENCES_AND_DISTINCT_GENE_RECORDS;FIRST_TWO_SELECTED_KEYS_NO_REORDER\n"));
		}
		System.err.println("NCRNA_SEED_SUPPORT_REPORTED");
	}

	private static LongList loadKeys(final String path, final int expected){
		final LongList keys=new LongList();
		final LongIntMap index=new LongIntMap(expected*2);
		final Throwable[] failure={null};
		RrnaResourceIO.stream(path, read->{
			if(failure[0]!=null){return;}
			try{
				require(read.bases!=null && read.bases.length==17, "Each ordered seed record must contain one complete 17-mer");
				long word=0;
				for(byte base:read.bases){
					final int code=AminoAcid.baseToNumber[base&255];
					require(code>=0, "Ordered seed keys cannot contain ambiguous bases");
					word=(word<<2)|code;
				}
				require(!index.contains(word), "Ordered seed keys must be unique");
				index.put(word, keys.size);
				keys.add(word);
			}catch(RuntimeException | Error error){
				failure[0]=error;
			}
		});
		rethrow(failure[0]);
		require(keys.size==expected, "The declared ordered-key count must match the complete file");
		return keys;
	}

	private static String word(final long key){
		final byte[] bases=new byte[17];
		for(int i=0; i<17; i++){bases[i]=(byte)"ACGT".charAt((int)((key>>>(2*(16-i)))&3));}
		return new String(bases, StandardCharsets.US_ASCII);
	}

	private static void rethrow(final Throwable error){
		if(error instanceof RuntimeException){throw (RuntimeException)error;}
		if(error instanceof Error){throw (Error)error;}
	}

	private static void require(final boolean value, final String why){
		if(!value){throw new IllegalArgumentException(why);}
	}

	private String seedPath;
	private String initialPath;
	private String out;
	private String[] genes;
	private int[] geneCounts;
	private int expectedInitial=-1;
	private int expectedKeys=-1;
	private final ByteBuilder row=new ByteBuilder();
	private static final long MASK=(1L<<34)-1;
}
