package test.ddl;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import ddl.DDLFamilyProfiles;
import ddl.DDLPentamerBlacklist;
import ddl.FamilySketchBuilder;
import ddl.ProteinKmerGeometry;
import dna.AminoAcid;

/** Protein-kmer geometry tests covering fractional and multiway reduced residues.
 * @author Yelan, Keqing
 */
public final class ProteinKmerGeometryTest{

	public static void main(String[] args){
		if(args.length!=0){throw new IllegalArgumentException("No arguments");}
		exactAndBinary();
		multiway();
		fullWordReduction();
		fullWordPositions();
		fullWordSketches();
		rejections();
		System.out.println("PROTEIN_KMER_GEOMETRY_TEST_PASS checks="+checks);
	}

	private static void exactAndBinary(){
		final long acdefg=rolling("ACDEFG");
		final long acdef=DDLPentamerBlacklist.pack("ACDEF".getBytes(StandardCharsets.US_ASCII));
		check(ProteinKmerGeometry.exact(5).word(acdef)==acdef, "Exact K5 preserves the raw five-bit protein word");
		final ProteinKmerGeometry k51=ProteinKmerGeometry.make(5, "hydro");
		check(k51.span()==6 && k51.suffixBits()==1, "K5.1 consumes one binary suffix residue");
		check(k51.word(acdefg)==((acdef<<1)|1), "Hydrophobic K5.1 appends the sixth-residue binary bit");
		final ProteinKmerGeometry k42=ProteinKmerGeometry.make(4, "hydro,hydro");
		check(k42.span()==6 && k42.suffixBits()==2, "K4.2 with comma partitions uses two suffix positions");
		final long acde=(acdef>>>5);
		check(k42.word(acdefg)==((acde<<2)|3), "K4.2 appends one bit for each of residues five and six");
		check(k42.pentamerPrefix(acdefg)==acdef, "K4.2 exposes the raw first pentamer for the existing blacklist");
		check(ProteinKmerGeometry.fromName(k42.name()).equals(k42), "Comma binary partition names round-trip");
	}

	private static void multiway(){
		final ProteinKmerGeometry explicit=ProteinKmerGeometry.make(4, "AGPST/RDQEK/NHFWY/CILMV");
		final ProteinKmerGeometry alias=ProteinKmerGeometry.make(4, "blosum4");
		check(alias.equals(explicit), "blosum4 is the configured four-way BLOSUM62 split");
		check(alias.span()==5 && alias.suffixBits()==2, "A four-way partition emits two bits at one suffix position");
		check((alias.word(rolling("AAAAA"))&3)==0, "BLOSUM group zero maps to bit pattern 0");
		check((alias.word(rolling("AAAAR"))&3)==1, "BLOSUM group one maps to bit pattern 1");
		check((alias.word(rolling("AAAAN"))&3)==2, "BLOSUM group two maps to bit pattern 2");
		check((alias.word(rolling("AAAAC"))&3)==3, "BLOSUM group three maps to bit pattern 3");
		check(ProteinKmerGeometry.fromName(alias.name()).equals(alias), "Slash multiway partition names round-trip");
		final ProteinKmerGeometry k44=ProteinKmerGeometry.make(4, "blosum4,blosum4");
		check(k44.span()==6 && k44.suffixBits()==4, "K4.4 uses two four-way suffix positions");
		check((k44.word(rolling("AAAARC"))&15)==7, "K4.4 packs each four-way suffix into two bits");
		check(ProteinKmerGeometry.make(4, "blosum").name().equals("exact4+ACILMFTWYV"), "blosum names the balanced two-way BLOSUM62 split");
	}

	private static void rejections(){
		reject(new Runnable(){@Override public void run(){ProteinKmerGeometry.make(4, "ACITV/LMFWY/RQEHK");}},
			"Multiway partitions must use a power-of-two group count");
		reject(new Runnable(){@Override public void run(){ProteinKmerGeometry.make(4, "ACITV/LMFWY/RQEHK/NDGP");}},
			"Multiway partitions must cover every canonical residue");
		reject(new Runnable(){@Override public void run(){ProteinKmerGeometry.make(4, "ACITV/LMFWY/RQEHKA/NDGPS");}},
			"Multiway partitions must not duplicate residues");
	}

	/** Literal symbol codes are independent of the optimized pair-lookup implementation. */
	private static void fullWordReduction(){
		final String groups="A/C/DE/FWY/G/IV/LM/NH/P/Q/RK/ST";
		final ProteinKmerGeometry k5=ProteinKmerGeometry.reduced(5,groups);
		final ProteinKmerGeometry k6=ProteinKmerGeometry.make(6,"all:"+groups);
		check(k5.k()==5 && k5.span()==5 && k5.exact()==0 && !k5.isExact(),
			"Fully reduced K5 consumes five residues without claiming an exact prefix");
		check(k6.k()==6 && k6.span()==6 && k6.suffixBits()==0 && !k6.isExact(6),
			"Fully reduced K6 is not a suffix-only or exact-six word");
		check(k5.word(rolling("ACDGF"))==0x01243L,
			"A,C,D,G,F have literal group codes0,1,2,4,3 in the frozen fixture partition");
		check(k6.word(rolling("ACDGFR"))==0x01243aL,
			"Reduced K6 packs all six four-bit symbols in sequence order");
		check(k6.word(rolling("ACEGYK"))==0x01243aL,
			"Substitutions within DE,FWY,RK merge at internal as well as terminal positions");
		check(k5.word(rolling("SSSSS"))==0xbbbbbL,
			"The twelfth symbol requires four bits; non-power-of-two alphabets must not truncate it");
		check(k6.rawPrefix(rolling("ACDGFR"),6)==rolling("ACDGFR"),
			"K6 raw prefix can expose the complete six-residue literal blacklist key");
		check(k6.pentamerPrefix(rolling("ACDGFR"))==rolling("ACDGF"),
			"Blacklist lookup retains the original raw pentamer, never the reduced word");
		check(k6.rollingMask()==((1L<<30)-1),
			"Raw six-residue rolling state still uses canonical five-bit residues");
		check(ProteinKmerGeometry.fromName(k5.name()).equals(k5) &&
			ProteinKmerGeometry.fromName(k6.name()).equals(k6),
			"Stored reduced geometry includes the complete alphabet and base word length");
		check(!k5.equals(k6) && !k5.equals(ProteinKmerGeometry.exact(5)),
			"Profile/query identity must distinguish alphabet and word length");
		reject(new Runnable(){@Override public void run(){k5.exactPrefix(rolling("ACDGF"));}},
			"A full-word reduction must not expose a nonexistent exact prefix");
		reject(new Runnable(){@Override public void run(){ProteinKmerGeometry.make(5,groups);}},
			"The legacy suffix grammar still rejects twelve groups without the explicit all prefix");
		reject(new Runnable(){@Override public void run(){ProteinKmerGeometry.reduced(6,"A/C/DE/FWY/G/IV/LM/NH/P/Q/RK");}},
			"Full-word reduction must cover all twenty standard residues");
		reject(new Runnable(){@Override public void run(){ProteinKmerGeometry.reduced(6,groups+"/A");}},
			"Full-word reduction must not assign a residue to two groups");
		reject(new Runnable(){@Override public void run(){ProteinKmerGeometry.reduced(13,groups);}},
			"Raw rolling state must fit the existing sixty-bit protein word contract");
		reject(new Runnable(){@Override public void run(){ProteinKmerGeometry.fromName("exact5+all:"+groups);}},
			"Stored exact-prefix identity must not be reinterpreted as a full-word reduction");
		reject(new Runnable(){@Override public void run(){ProteinKmerGeometry.fromName("exact5+ all:"+groups);}},
			"Whitespace must not bypass the stored exact-prefix identity check");
		reject(new Runnable(){@Override public void run(){ProteinKmerGeometry.fromName("reduced5");}},
			"Stored reduced geometry must bind its alphabet explicitly");
	}

	/** Every position uses the same literal twelve-symbol alphabet, including odd leading residues. */
	private static void fullWordPositions(){
		final String alphabet="ARNDCQEGHILKMFPSTWYV";
		final int[] expected={0,10,7,2,1,9,2,4,7,5,6,10,6,3,8,11,11,3,3,5};
		for(int k=5; k<=6; k++){
			final ProteinKmerGeometry geometry=ProteinKmerGeometry.reduced(k,"A/C/DE/FWY/G/IV/LM/NH/P/Q/RK/ST");
			final char[] word=new char[k];Arrays.fill(word,'A');
			for(int position=0; position<k; position++){
				for(int residue=0; residue<alphabet.length(); residue++){
					word[position]=alphabet.charAt(residue);
					final long literal=((long)expected[residue])<<(4*(k-1-position));
					check(geometry.word(rolling(new String(word)))==literal,
						"Every raw residue must map at its own position: k="+k+" position="+position+" residue="+word[position]);
				}
				word[position]='A';
			}
		}
	}

	/** Existing production sketch loops consume projected words without rewriting biological sequences. */
	private static void fullWordSketches(){
		final ProteinKmerGeometry geometry=ProteinKmerGeometry.reduced(6,"A/C/DE/FWY/G/IV/LM/NH/P/Q/RK/ST");
		final byte[] raw="ACDGFRXACEGYK".getBytes(StandardCharsets.US_ASCII), original=raw.clone();
		final byte[] encoded=new byte[raw.length];
		for(int i=0; i<raw.length; i++){
			encoded[i]=(byte)(raw[i]=='X' ? AminoAcid.acidToNumberExtended['X'] : AminoAcid.acidToNumberNoStops[raw[i]]);
		}
		final char[] direct=new char[512], encodedSketch=new char[512], profile=new char[512];
		check(DDLFamilyProfiles.sketch(raw,geometry,0,direct,null)==2 &&
			DDLFamilyProfiles.sketchEncoded(encoded,geometry,0,encodedSketch,null)==2,
			"Ambiguity splits two six-residue words and contributes no cross-boundary word");
		check(Arrays.equals(direct,encodedSketch),"ASCII and encoded reduced-word sketches must agree");
		int occupied=0;for(char value:direct){if(value!=0){occupied++;}}
		check(occupied==1,"The two distinct raw words have the same literal reduced code and one plain register winner");
		final FamilySketchBuilder builder=new FamilySketchBuilder(512,geometry,0,false,1,0,null);
		builder.reset();builder.addRaw(raw);builder.writeProfile(profile);
		check(builder.members()==1 && builder.validWords()==2 && Arrays.equals(direct,profile),
			"A one-member modal profile must retain the member's complete reduced sketch");
		check(DDLFamilyProfiles.sketchPositional(raw,geometry,0,4,2,direct,null)==2 &&
			DDLFamilyProfiles.sketchPositionalEncoded(encoded,geometry,0,4,2,encodedSketch,null)==2 &&
			Arrays.equals(direct,encodedSketch),"Positional member/query paths share reduction and raw span boundaries");
		check(Arrays.equals(raw,original),"Sketch projection must not modify the biological amino-acid input");
	}

	private static long rolling(final String sequence){
		long word=0;
		for(int i=0; i<sequence.length(); i++){
			final int aa=AminoAcid.acidToNumberNoStops[sequence.charAt(i)];
			check(aa>=0 && aa<20, "Test sequences use only standard amino acids");
			word=(word<<5)|aa;
		}
		return word;
	}

	private static void reject(Runnable task, String message){
		checks++;
		try{task.run();}catch(IllegalArgumentException expected){return;}
		throw new AssertionError(message);
	}

	private static void check(boolean condition, String message){
		checks++;
		if(!condition){throw new AssertionError(message);}
	}

	private static long checks;
}
