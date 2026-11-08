package test.bloom;

import java.io.BufferedWriter;
import java.io.IOException;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;

import bloom.BloomFilter;
import bloom.KCountArray;
import bloom.KCountArray7MTA;
import stream.Read;
import structures.IntList;
import structures.LongList;
import ukmer.Kmer;

/** Finite query/scratch-buffer contracts and explicit diagnostics for ambiguous policies.
 * @author Brian Bushnell, Shinobu, Jean, Sucrose */
public final class BloomFilterContractTest{

	/** Runs one named case in its own process so every baseline outcome is retained. */
	public static void main(final String[] args) throws Exception{
		String name=null;
		Path outdir=null;
		for(String arg : args){
			if(arg.startsWith("case=")){name=arg.substring(5);}
			else if(arg.startsWith("outdir=")){outdir=Paths.get(arg.substring(7));}
			else{throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		if(name==null){
			throw new IllegalArgumentException("Specify case=streak|short|big-short|big-undefined|"+
					"average-diagnostic|fraction-diagnostic|constructor-primary|constructor-primary-equals|"+
					"constructor-extra|constructor-extra-equals|fraction-short|fraction-long|fraction-continuation|"+
					"average-endpoints|rounding-arithmetic");
		}
		final boolean constructorCase=name.startsWith("constructor-");
		final boolean fractionCase=name.equals("fraction-short") || name.equals("fraction-long") || name.equals("fraction-continuation");
		final boolean averageCase=name.equals("average-endpoints");
		final boolean arithmeticCase=name.equals("rounding-arithmetic");
		if((constructorCase || fractionCase || averageCase || arithmeticCase)!=(outdir!=null)){
			throw new IllegalArgumentException("Constructor/fraction/endpoint/arithmetic cases require outdir=<new directory>.");
		}
		final long oldCells=BloomFilter.OVERRIDE_CELLS;
		final boolean oldPrint=BloomFilter.printMem;
		try{
			BloomFilter.OVERRIDE_CELLS=100003;
			BloomFilter.printMem=false;
			KCountArray7MTA.setSeed(917);
			if(name.equals("streak")){streak();}
			else if(name.equals("short")){shortCounts();}
			else if(name.equals("big-short")){shortBigCounts();}
			else if(name.equals("big-undefined")){undefinedBigCounts();}
			else if(name.equals("average-diagnostic")){averageDiagnostic();}
			else if(name.equals("fraction-diagnostic")){fractionDiagnostic();}
			else if(name.equals("fraction-short")){definedFraction(outdir, 31);}
			else if(name.equals("fraction-long")){definedFraction(outdir, 63);}
			else if(name.equals("fraction-continuation")){definedFraction(outdir, 31, true);}
			else if(name.equals("average-endpoints")){averageEndpoints(outdir);}
			else if(name.equals("rounding-arithmetic")){roundingArithmetic(outdir);}
			else if(name.equals("constructor-primary")){constructorPath(name, outdir, true, false);}
			else if(name.equals("constructor-primary-equals")){constructorPath(name, outdir, true, true);}
			else if(name.equals("constructor-extra")){constructorPath(name, outdir, false, false);}
			else if(name.equals("constructor-extra-equals")){constructorPath(name, outdir, false, true);}
			else{throw new IllegalArgumentException("Unknown case: "+name);}
			System.out.println((name.endsWith("-diagnostic") ?
					"BLOOM_FILTER_DIAGNOSTIC_COMPLETE case=" : "BLOOM_FILTER_CONTRACT_PASS case=")+name);
		}finally{
			BloomFilter.OVERRIDE_CELLS=oldCells;
			BloomFilter.printMem=oldPrint;
		}
	}

	/** A run cannot span an undefined base; both public query overloads must agree. */
	private static void streak(){
		final BloomFilter filter=filter(3, false);
		filter.filter.increment(6, 1);//ACG
		filter.filter.increment(27, 1);//CGT
		final Read contiguous=read("ACGT");
		check(!filter.passes(contiguous, 1), "Two adjacent supported windows must match.");
		check(!filter.passes(contiguous, new LongList(), 1), "Buffered control must match.");
		final Read separated=read("ACGNACG");
		final boolean direct=filter.passes(separated, 1);
		final boolean buffered=filter.passes(separated, new LongList(), 1);
		System.out.println("separated directPass="+direct+" bufferedPass="+buffered);
		check(direct && buffered, "Two isolated hits separated by N are not consecutive windows.");
	}

	/** Reusing a count list for a short sequence must not retain the preceding sequence. */
	private static void shortCounts(){
		final BloomFilter filter=filter(3, false);
		filter.filter.increment(6, 1);
		final IntList counts=new IntList();
		check(filter.fillCounts(bytes("ACGN"), counts)==1 && counts.size==2,
				"A full input must emit one count per window, including zero for undefined windows.");
		check(counts.get(0)==1 && counts.get(1)==0, "Control counts differ from the explicit sequence.");
		final int valid=filter.fillCounts(bytes("AC"), counts);
		System.out.println("short valid="+valid+" size="+counts.size);
		check(valid==0 && counts.size==0, "A sequence shorter than k has no output windows.");
	}

	/** Sub-kmer counts are not full-sized windows when the requested big k does not fit. */
	private static void shortBigCounts(){
		final BloomFilter filter=filter(5, false);
		filter.filter.increment(6, 1);
		filter.filter.increment(27, 1);
		final IntList counts=new IntList();
		counts.add(99);
		final int valid=filter.fillCountsBig(bytes("ACGT"), counts);
		System.out.println("big-short valid="+valid+" size="+counts.size);
		check(valid==0 && counts.size==0, "Four bases cannot produce a five-base window.");
	}

	/** Checks full-window geometry by independently examining every base in each window. */
	private static void undefinedBigCounts(){
		final BloomFilter filter=filter(5, false);
		for(int key=0; key<64; key++){filter.filter.increment(key, 1);}
		for(int key=0; key<64; key++){
			check(filter.getCount(key)==1, "A hash collision would invalidate the uniform count-one fixture.");
		}
		int failures=0;
		for(String sequence : new String[]{"ACGTNACGT", "NNNNNNNN", "ACGTACGT", "ACGTNACGTAC"}){
			final byte[] bases=bytes(sequence);
			final int windows=bases.length-4;
			final int[] expected=new int[windows];
			int expectedValid=0;
			for(int start=0; start<windows; start++){
				boolean defined=true;
				for(int pos=start; pos<start+5; pos++){if(bases[pos]=='N'){defined=false;}}
				expected[start]=(defined ? 1 : 0);
				if(defined){expectedValid++;}
			}
			final IntList counts=new IntList();
			counts.add(99);
			final int valid=filter.fillCountsBig(bases, counts);
			boolean pass=valid==expectedValid && counts.size==windows;
			for(int i=0; i<windows; i++){pass&=counts.get(i)==expected[i];}
			System.out.println("sequence="+sequence+" expectedValid="+expectedValid+" valid="+valid+
					" expectedSize="+windows+" size="+counts.size+" counts="+counts+" pass="+pass);
			if(!pass){failures++;}
		}
		check(failures==0, "Big-window counts and validity must respect every undefined base; failures="+failures);
	}

	/** Reports edge-smoothing orientation behavior without choosing a replacement policy. */
	private static void averageDiagnostic(){
		final BloomFilter filter=filter(3, true);
		filter.filter.increment(63, 10);//Canonical AAA/TTT
		filter.filter.increment(47, 1);//Canonical AAC/GTT
		check(filter.getCount(63)==10 && filter.getCount(47)==1,
				"The finite table must retain distinct seeded control counts.");
		final float forward=filter.averageCount(bytes("AAAC"));
		final float reverse=filter.averageCount(bytes("GTTT"));
		System.out.println("average forward="+forward+" reverseComplement="+reverse);
	}

	/** Reports undefined-window handling without deciding the fraction denominator policy. */
	private static void fractionDiagnostic(){
		final BloomFilter filter=filter(3, false);
		final Read undefined=read("NNNN");
		final boolean passes=filter.hasHighCountFraction(undefined, 1, 1f);
		final float fraction=filter.highCountFraction(undefined, 1, false);
		System.out.println("undefined hasHighFraction="+passes+" numericHighFraction="+fraction);
	}

	/** Exercises the defined-window policy in the actual short or multiword query path. */
	private static void definedFraction(final Path directory, final int k) throws IOException{
		definedFraction(directory, k, false);
	}

	/** Shares verified table setup between the original policy rows and separate continuation controls. */
	private static void definedFraction(final Path directory, final int k, final boolean continuation) throws IOException{
		Files.createDirectory(directory);
		final boolean oldPacked=Kmer.PACKED, oldMix=Kmer.FULL_MIX, oldMask=Kmer.MASK_CORE;
		try{
			Kmer.PACKED=Kmer.FULL_MIX=true;
			Kmer.MASK_CORE=false;
			check(Kmer.getKbig(k)==k, "The fraction fixture must represent the requested full kmer span.");
			final BloomFilter filter=new BloomFilter(k, k, 8, 1, 1, true, false, false, 0.01f);
			final String high=repeat('A', k), low=repeat('C', k);
			final long highKey=fractionKey(high, k), lowKey=fractionKey(low, k);
			final long highCell=physicalLane(filter.filter, highKey), lowCell=physicalLane(filter.filter, lowKey);
			filter.filter.increment(highKey, 2);
			filter.filter.increment(lowKey, 1);
			check(highCell!=lowCell && filter.getCount(highKey)==2 && filter.getCount(lowKey)==1,
					"Depth2 high and depth1 low seeds must occupy separate cells in this table.");
			//Arithmetic control only: no allocation or execution of a 16-million-window read.
			final int possible=16777222, defined=16777221;
			check(possible-Math.round(1f*possible)==0 && defined-Math.round(1f*defined)==1,
					"Binary32 rounding above 2^24 makes the possible-window budget nonmonotone.");
			int failures=0;
			try(BufferedWriter writer=Files.newBufferedWriter(directory.resolve("fraction.tsv"), StandardCharsets.US_ASCII)){
				writer.write("SETUP k="+k+" represented="+Kmer.getKbig(k)+" words="+Kmer.getMult(k)+
						" cells="+filter.filter.cells+" highKey="+highKey+" lowKey="+lowKey+
						" highCell="+highCell+" lowCell="+lowCell+" highDepth=2 lowDepth=1 threshold=2"+
						" packed=true fullMix=true maskCore=false arithmeticPossible=16777222"+
						" arithmeticDefined=16777221 possibleBudget=0 definedBudget=1\n");
				writer.write("case\tk\tbases\tfraction\tdefined\thigh\troundedRequired\texpected\tactual\tinputPreserved\n");
				if(continuation){
					failures+=fractionRow(writer, filter, k, "trailing-N", high+'N', 1f, 1, 1, true);
					failures+=fractionRow(writer, filter, k, "leading-N", 'N'+high, 1f, 1, 1, true);
					failures+=fractionRow(writer, filter, k, "consecutive-N", high+"NN"+low, .75f, 2, 1, false);
					final String highLong=high+'A', lowLong=low+'C', mixed=highLong+'N'+lowLong;
					failures+=fractionRow(writer, filter, k, "long-prefix-half", mixed, .5f, 4, 2, true);
					failures+=fractionRow(writer, filter, k, "long-prefix-rounding", mixed, .75f, 4, 2, false);
					failures+=fractionRow(writer, filter, k, "carried-low-two", lowLong+'N'+highLong, .75f, 4, 2, false);
				}else{
					failures+=fractionRow(writer, filter, k, "all-high", high+'A', 1f, 2, 2, true);
					failures+=fractionRow(writer, filter, k, "all-low", low+'C', .5f, 2, 0, false);
					failures+=fractionRow(writer, filter, k, "zero-fraction", low, 0f, 1, 0, true);
					final String mixed=high+'N'+low;
					failures+=fractionRow(writer, filter, k, "mixed-half", mixed, .5f, 2, 1, true);
					failures+=fractionRow(writer, filter, k, "mixed-below-rounding", mixed, Math.nextDown(.75f), 2, 1, true);
					failures+=fractionRow(writer, filter, k, "mixed-rounding", mixed, .75f, 2, 1, false);
					failures+=fractionRow(writer, filter, k, "mixed-one", mixed, 1f, 2, 1, false);
					failures+=fractionRow(writer, filter, k, "mixed-reversed", low+'N'+high, .75f, 2, 1, false);
					final String undefined=repeat('N', k+1), shortRead=repeat('A', k-1);
					failures+=fractionRow(writer, filter, k, "all-N-zero", undefined, 0f, 0, 0, false);
					failures+=fractionRow(writer, filter, k, "all-N-one", undefined, 1f, 0, 0, false);
					failures+=fractionRow(writer, filter, k, "short-zero", shortRead, 0f, 0, 0, false);
					failures+=fractionRow(writer, filter, k, "short-one", shortRead, 1f, 0, 0, false);
					failures+=fractionRow(writer, filter, k, "shorter", repeat('A', k-2), 1f, 0, 0, false);
					failures+=fractionRow(writer, filter, k, "short-islands", shortRead+'N'+repeat('C', k-1), 1f, 0, 0, false);
					failures+=fractionRow(writer, filter, k, "null", null, 1f, 0, 0, false);
				}
			}
			check(failures==0, "Defined-window fraction failed "+failures+" of "+(continuation ? 6 : 15)+" rows at k="+k);
		}finally{
			Kmer.PACKED=oldPacked; Kmer.FULL_MIX=oldMix; Kmer.MASK_CORE=oldMask;
		}
	}

	/** Checks literal window totals independently, then records the actual predicate result. */
	private static int fractionRow(final BufferedWriter writer, final BloomFilter filter, final int k,
			final String name, final String sequence, final float fraction, final int expectedDefined,
			final int expectedHigh, final boolean expected) throws IOException{
		int defined=0, high=0;
		for(int start=0; sequence!=null && start+k<=sequence.length(); start++){
			boolean allA=true, allC=true, valid=true;
			for(int i=start; i<start+k; i++){
				final char base=sequence.charAt(i);
				allA&=base=='A'; allC&=base=='C'; valid&=base=='A' || base=='C';
			}
			if(valid){
				check(allA || allC, "Fixture window must be one of the two independently seeded homopolymers.");
				defined++;
				if(allA){high++;}
			}
		}
		check(defined==expectedDefined && high==expectedHigh,
				"Literal fixture window totals changed for "+name+": defined="+defined+" high="+high);
		final Read query=(sequence==null ? null : read(sequence));
		final boolean actual=filter.hasHighCountFraction(query, 2, fraction);
		final boolean preserved=query==null || Arrays.equals(bytes(sequence), query.bases);
		final String row=name+'\t'+k+'\t'+(sequence==null ? "<null>" : sequence)+'\t'+fraction+'\t'+
				defined+'\t'+high+'\t'+Math.round(fraction*defined)+'\t'+expected+'\t'+actual+'\t'+preserved;
		writer.write(row); writer.newLine();
		System.out.println(row);
		return actual==expected && preserved ? 0 : 1;
	}

	/** Encodes one known homopolymer with the matching short or full-span public key API. */
	private static long fractionKey(final String sequence, final int k){
		assert(sequence.length()==k) : "Seed must describe exactly one full window; length="+sequence.length()+" k="+k;
		if(k<=31){
			long key=0;
			for(int i=0; i<k; i++){key=(key<<2)|(sequence.charAt(i)=='A' ? 0 : 1);}
			return KCountArray.makeCanonical2(key, k);
		}
		final Kmer key=new Kmer(k);
		for(int i=0; i<k; i++){key.addRight(sequence.charAt(i));}
		return key.xor();
	}

	/** Builds a short literal homopolymer for independent finite window counts. */
	private static String repeat(final char base, final int length){
		final char[] bases=new char[length];
		Arrays.fill(bases, base);
		return new String(bases);
	}

	/** Checks both endpoint orientations while retaining singleton, interior and flat controls. */
	private static void averageEndpoints(final Path directory) throws IOException{
		Files.createDirectory(directory);
		final BloomFilter filter=filter(3, true);
		final long[] keys={63, 47, 62};//AAA/TTT, AAC/GTT, CAA/TTG.
		final int[] depths={10, 1, 1};
		final long[] cells=new long[keys.length];
		for(int i=0; i<keys.length; i++){
			cells[i]=physicalLane(filter.filter, keys[i]);
			for(int j=0; j<i; j++){check(cells[i]!=cells[j], "Endpoint fixture keys must occupy distinct cells.");}
			filter.filter.increment(keys[i], depths[i]);
		}
		for(int i=0; i<keys.length; i++){
			check(filter.getCount(keys[i])==depths[i], "Endpoint seed count changed for key "+keys[i]);
		}
		int failures=0;
		try(BufferedWriter writer=Files.newBufferedWriter(directory.resolve("average.tsv"), StandardCharsets.US_ASCII)){
			writer.write("SETUP k=3 canonical=true keys="+Arrays.toString(keys)+" depths="+Arrays.toString(depths)+
					" cells="+Arrays.toString(cells)+"\n");
			writer.write("case\tbases\tdepths\texpected\tactual\tinputPreserved\n");
			failures+=averageRow(writer, filter, "left-spike", "AAAC", new int[]{10, 1}, 1f);
			failures+=averageRow(writer, filter, "right-spike-reverse", "GTTT", new int[]{1, 10}, 1f);
			failures+=averageRow(writer, filter, "singleton", "AAA", new int[]{10}, 10f);
			failures+=averageRow(writer, filter, "interior", "CAAAC", new int[]{1, 10, 1}, 1f);
			failures+=averageRow(writer, filter, "interior-reverse", "GTTTG", new int[]{1, 10, 1}, 1f);
			failures+=averageRow(writer, filter, "flat", "AAAA", new int[]{10, 10}, 10f);
			failures+=averageRow(writer, filter, "all-N", "NNNN", new int[0], 0f);
			failures+=averageRow(writer, filter, "short", "AA", new int[0], 0f);
			failures+=averageRow(writer, filter, "null", null, new int[0], 0f);
		}
		check(failures==0, "Endpoint smoothing failed "+failures+" of 9 rows.");
	}

	/** Enumerates literal three-base keys independently and records the smoothing result. */
	private static int averageRow(final BufferedWriter writer, final BloomFilter filter,
			final String name, final String sequence, final int[] depths, final float expected) throws IOException{
		int counted=0;
		for(int start=0; sequence!=null && start+3<=sequence.length(); start++){
			long key=0, reverse=0;
			boolean defined=true;
			for(int offset=0; offset<3; offset++){
				final int base="ACGT".indexOf(sequence.charAt(start+offset));
				if(base<0){defined=false; break;}
				key=(key<<2)|base;
				reverse|=(long)(3-base)<<(2*offset);
			}
			if(defined){
				check(counted<depths.length && filter.getCount(Math.max(key, reverse))==depths[counted],
						"Literal endpoint depth sequence changed for "+name+" at window "+counted);
				counted++;
			}
		}
		check(counted==depths.length, "Defined endpoint window count changed for "+name);
		final byte[] bases=(sequence==null ? null : bytes(sequence));
		final float actual=filter.averageCount(bases);
		final boolean preserved=bases==null || Arrays.equals(bytes(sequence), bases);
		final String row=name+'\t'+(sequence==null ? "<null>" : sequence)+'\t'+Arrays.toString(depths)+
				'\t'+expected+'\t'+actual+'\t'+preserved;
		writer.write(row); writer.newLine();
		System.out.println(row);
		return actual==expected && preserved ? 0 : 1;
	}

	/** Compares the actual private rounding helper against Math.round on a finite binary32 grid.
	 * Large window counts are arithmetic arguments only; no corresponding reads are allocated. */
	private static void roundingArithmetic(final Path directory) throws Exception{
		Files.createDirectory(directory);
		final Method helper=BloomFilter.class.getDeclaredMethod("requiredHighCount", float.class, int.class);
		helper.setAccessible(true);
		final int[] windows={0, 1, 2, 3, 31, 127, 128, 150, 1000,
				(1<<23)-1, 1<<23, (1<<23)+1, (1<<24)-1, 1<<24, (1<<24)+1, Integer.MAX_VALUE};
		final float[] fractions={Float.NEGATIVE_INFINITY, -Float.MAX_VALUE, -2f, -1.5f, -1f,
				Math.nextDown(-.5f), -.5f, Math.nextUp(-.5f), -Float.MIN_NORMAL,
				-Math.nextDown(Float.MIN_NORMAL), -Float.MIN_VALUE, -0f, 0f, Float.MIN_VALUE,
				Math.nextDown(Float.MIN_NORMAL), Float.MIN_NORMAL,
				Math.nextDown(.5f), .5f, Math.nextUp(.5f), Math.nextDown(.75f), .75f, Math.nextUp(.75f),
				Math.nextDown(1f), 1f, Math.nextUp(1f), 1.5f, 2f, Float.MAX_VALUE, Float.POSITIVE_INFINITY, Float.NaN};
		int failures=0, rows=0;
		try(BufferedWriter writer=Files.newBufferedWriter(directory.resolve("rounding.tsv"), StandardCharsets.US_ASCII)){
			writer.write("fractionBits\tfraction\twindows\tproductBits\tproduct\texpected\tactual\terror\tpass\n");
			for(final float fraction : fractions){
				for(final int count : windows){
					final float product=fraction*count;
					final int expected=Math.round(fraction*count);
					int actual=0;
					String error="none";
					try{
						actual=((Integer)helper.invoke(null, fraction, count)).intValue();
					}catch(InvocationTargetException e){
						error=e.getCause().getClass().getName();
					}
					final boolean pass=error.equals("none") && actual==expected;
					writer.write(Integer.toHexString(Float.floatToRawIntBits(fraction))+"\t"+fraction+"\t"+count+"\t"+
							Integer.toHexString(Float.floatToRawIntBits(product))+"\t"+product+"\t"+
							expected+"\t"+actual+"\t"+error+"\t"+pass+"\n");
					rows++;
					if(!pass){failures++;}
				}
			}
		}
		System.out.println("ROUNDING_ARITHMETIC rows="+rows+" failures="+failures+" oracle=Math.round(floatProduct)");
		check(failures==0, "requiredHighCount disagreed with Math.round on "+failures+" of "+rows+" arithmetic rows.");
	}

	/** Equals-containing constructor paths must preserve everything after the first delimiter. */
	private static void constructorPath(final String name, final Path directory, final boolean primary,
			final boolean equalsPath) throws IOException{
		Files.createDirectory(directory);
		final String label=(primary ? "primary" : "extra");
		final Path decoy=directory.resolve(label+"-decoy.fq");
		final Path target=directory.resolve(equalsPath ? label+"-decoy.fq=chosen.fq" : label+"-chosen.fq");
		final Path primaryPath=primary ? target : directory.resolve("ordinary-primary.fq");
		final byte[] primaryBytes=fastq("primary", "ACG");
		final byte[] targetBytes=(primary ? primaryBytes : fastq("target", "GTA"));
		final byte[] decoyBytes=fastq("decoy", "CCC");
		final byte[] expectedPrimaryBytes=(primary ? targetBytes : primaryBytes);
		Files.write(target, primary ? primaryBytes : targetBytes);
		Files.write(decoy, decoyBytes);
		if(!primary){Files.write(primaryPath, primaryBytes);}
		final BloomFilter filter=new BloomFilter(primary ?
				new String[]{"extra=", "in="+target, "k=3", "kbig=3", "bits=8", "hashes=1", "rcomp=f"} :
				new String[]{"extra=", "in="+primaryPath, "extra="+target, "k=3", "kbig=3", "bits=8", "hashes=1", "rcomp=f"});
		final int actualAcg=filter.getCount(6);
		final int actualGta=filter.getCount(44);
		final int actualCcc=filter.getCount(21);
		final int expectedAcg=1;
		final int expectedGta=(primary ? 0 : 1);
		final int expectedCcc=0;
		final long cellAcg=physicalLane(filter.filter, 6);
		final long cellGta=physicalLane(filter.filter, 44);
		final long cellCcc=physicalLane(filter.filter, 21);
		final boolean cellsSeparate=cellAcg!=cellGta && cellAcg!=cellCcc && cellGta!=cellCcc;
		final boolean primaryPreserved=Arrays.equals(expectedPrimaryBytes, Files.readAllBytes(primaryPath));
		final boolean targetPreserved=Arrays.equals(targetBytes, Files.readAllBytes(target));
		final boolean decoyPreserved=Arrays.equals(decoyBytes, Files.readAllBytes(decoy));
		final String result="case\t"+name+"\tprimary\t"+primary+"\tequalsPath\t"+equalsPath+
				"\texpected6\t"+expectedAcg+"\tactual6\t"+actualAcg+
				"\texpected44\t"+expectedGta+"\tactual44\t"+actualGta+
				"\texpected21\t"+expectedCcc+"\tactual21\t"+actualCcc+
				"\tcell6\t"+cellAcg+"\tcell44\t"+cellGta+"\tcell21\t"+cellCcc+
				"\tcellsSeparate\t"+cellsSeparate+"\tprimaryPreserved\t"+primaryPreserved+
				"\ttargetPreserved\t"+targetPreserved+"\tdecoyPreserved\t"+decoyPreserved;
		System.out.println(result);
		Files.write(directory.resolve("observed.tsv"), bytes(result+"\n"));
		check(cellsSeparate, "Fixture keys 6/44/21 must not alias in the loaded table.");
		check(primaryPreserved, "Primary input bytes changed.");
		check(targetPreserved, "Target input bytes changed.");
		check(decoyPreserved, "Decoy input bytes changed.");
		check(actualAcg==expectedAcg, "ACG count must distinguish the ordinary primary from the decoy.");
		check(actualGta==expectedGta, "GTA count must distinguish the extra target from the decoy.");
		check(actualCcc==expectedCcc, "The valid decoy prefix file must not be loaded.");
	}

	/** Uses a tiny fixed-capacity, one-hash table; ksmall is always three. */
	private static BloomFilter filter(final int bigK, final boolean canonical){
		return new BloomFilter(3, bigK, 8, 1, 2, canonical, false, false, 0.01f);
	}

	/** Returns the row-major physical cell for a one-lane raw key in this table. */
	private static long physicalLane(final KCountArray7MTA array, final long rawKey) throws IOException{
		long key=array.hashForLane(rawKey, 0, 0);
		final long arrayMask=longField(ARRAY_MASK, array);
		final int arrayBits=(int)longField(ARRAY_BITS, array);
		final long cellMod=longField(CELL_MOD, array);
		final long num=key&arrayMask;
		key=(key>>>arrayBits)%cellMod;
		return num*cellMod+key;
	}

	/** Makes an unpaired, unmodified sequence for query-only checks. */
	private static Read read(final String bases){
		return new Read(bytes(bases), null, "fixture", 0, false);
	}

	/** Encodes the literal ASCII test sequences. */
	private static byte[] bytes(final String value){return value.getBytes(StandardCharsets.US_ASCII);}

	/** Encodes a single-read FASTQ with the supplied bases. */
	private static byte[] fastq(final String name, final String bases){
		return bytes("@"+name+"\n"+bases+"\n+\nIII\n");
	}

	/** Keeps failed contracts fatal even with JVM assertions disabled. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Reads an integer-like reflected field as a long. */
	private static long longField(final Field field, final Object owner) throws IOException{
		try{
			return ((Number)field.get(owner)).longValue();
		}catch(IllegalAccessException e){
			throw new IOException(e);
		}
	}

	private static final Field ARRAY_MASK;
	private static final Field ARRAY_BITS;
	private static final Field CELL_MOD;

	static{
		try{
			ARRAY_MASK=bloom.KCountArray.class.getDeclaredField("arrayMask");
			ARRAY_BITS=bloom.KCountArray.class.getDeclaredField("arrayBits");
			CELL_MOD=KCountArray7MTA.class.getDeclaredField("cellMod");
			ARRAY_MASK.setAccessible(true);
			ARRAY_BITS.setAccessible(true);
			CELL_MOD.setAccessible(true);
		}catch(Exception e){
			throw new RuntimeException(e);
		}
	}
}
