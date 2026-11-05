package test.ukmer;

import java.lang.reflect.Field;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.Arrays;

import kmer.AbstractKmerTableSet;
import ukmer.Kmer;
import ukmer.KmerTableSetU;

/**
 * Finite long-kmer extra-file expansion checks through construction and real
 * input loading.  Run one case per fresh JVM, with a separate dir for each
 * assertion mode/case.  Existing input files are never overwritten.
 * @author Jean
 */
public final class KmerExtraInputLongTest {

	/** Accepts case=expanded|explicit|literal-hash and dir= for retained inputs. */
	public static void main(String[] args) throws Exception{
		String mode=null;
		Path dir=null;
		for(String arg : args){
			if(arg.startsWith("case=")){mode=arg.substring(5);}
			else if(arg.startsWith("dir=")){dir=Paths.get(arg.substring(4));}
			else{throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		if(dir==null || !("expanded".equals(mode) || "explicit".equals(mode) || "literal-hash".equals(mode))){
			throw new IllegalArgumentException("Expected case=expanded|explicit|literal-hash and dir=<new input directory>");
		}

		Files.createDirectories(dir);
		final String mainBases=repeat("ACGTAA", 10)+"ACG";
		final String firstBases=repeat("CCGAATG", 9);
		final String secondBases=repeat("GATACCA", 9);
		final Path main=dir.resolve("main.fq"), first=dir.resolve("extra1.fq"), second=dir.resolve("extra2.fq");
		write(main, "main", mainBases);
		write(first, "extra1", firstBases);
		write(second, "extra2", secondBases);

		final boolean literal=mode.equals("literal-hash");
		final String extra;
		final String[] expectedBases;
		if(literal){
			final Path path=dir.resolve("literal#.fq");
			write(path, "literal", firstBases);
			extra=path.toString();
			expectedBases=new String[] {mainBases, firstBases};
		}else if(mode.equals("expanded")){
			extra=dir.resolve("extra#.fq").toString();
			expectedBases=new String[] {mainBases, firstBases, secondBases};
		}else{
			extra=first+","+second;
			expectedBases=new String[] {mainBases, firstBases, secondBases};
		}

		final Field rcomp=Kmer.class.getDeclaredField("rcomp");
		rcomp.setAccessible(true);
		final boolean oldPacked=Kmer.PACKED, oldRcomp=rcomp.getBoolean(null);
		final boolean oldMaskCore=AbstractKmerTableSet.MASK_CORE;
		final boolean oldMaskMiddle=AbstractKmerTableSet.MASK_MIDDLE;
		try{
			Kmer.PACKED=true;
			AbstractKmerTableSet.MASK_CORE=false;
			AbstractKmerTableSet.MASK_MIDDLE=false;
			final KmerTableSetU tables=new KmerTableSetU(new String[] {
				"in="+main, "extra="+extra, "k=63", "ways=1", "threads=1", "prealloc=f",
				"prefilter=0", "onepass=f", "rcomp=t", "qin=33", "interleaved=f",
				"minprob=0", "minavgquality=0", "qtrim=f", "ecco=f", "merge=f"
			}, 12);
			final long added=tables.processInput();
			check(!tables.errorState, "Input loading reported an error");
			check(tables.kbig()==63 && tables.readsIn==expectedBases.length && tables.basesIn==63L*expectedBases.length,
				"Wrong long-kmer geometry or input read/base totals");
			check(added==expectedBases.length && tables.kmersLoaded==expectedBases.length &&
					tables.kmersIn==expectedBases.length,
				"Every 63-base read must contribute one distinct long k-mer");
			for(int i=0; i<expectedBases.length; i++){
				for(int j=i+1; j<expectedBases.length; j++){
					check(!Arrays.equals(kmer(expectedBases[i]).key(), kmer(expectedBases[j]).key()),
						"Fixture query keys must not collapse.");
				}
			}
			for(String bases : expectedBases){
				check(tables.getCount(kmer(bases))==1, "Wrong retained count for "+bases);
			}
			check(tables.walk().next(), "Loaded table must be walkable after real input processing");
		}finally{
			AbstractKmerTableSet.MASK_MIDDLE=oldMaskMiddle;
			AbstractKmerTableSet.MASK_CORE=oldMaskCore;
			rcomp.setBoolean(null, oldRcomp);
			Kmer.PACKED=oldPacked;
		}
		System.out.println("PASS case="+mode+" reads="+expectedBases.length+" distinctKmers="+expectedBases.length);
	}

	/** Builds a 63-base packed query through the ordinary rolling API. */
	private static Kmer kmer(final String bases){
		final Kmer kmer=new Kmer(bases.length());
		for(int i=0; i<bases.length(); i++){kmer.addRight(bases.charAt(i));}
		return kmer;
	}

	/** Repeats a short motif into one retained test sequence. */
	private static String repeat(final String motif, final int count){
		final StringBuilder sb=new StringBuilder(motif.length()*count);
		for(int i=0; i<count; i++){sb.append(motif);}
		return sb.toString();
	}

	/** Writes one high-quality unpaired read, retaining exact input bytes. */
	private static void write(final Path path, final String id, final String bases) throws Exception{
		Files.write(path, ("@"+id+"\n"+bases+"\n+\n"+repeat("I", bases.length())+"\n").
				getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
	}

	/** Keeps the oracle active with Java assertions disabled. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
}
