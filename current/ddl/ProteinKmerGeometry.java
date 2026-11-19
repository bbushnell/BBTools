package ddl;

import java.util.Arrays;
import java.util.Locale;

import dna.AminoAcid;

/** Protein words with exact prefixes/reduced suffixes, full-word alphabets or ordered dimer classes.
 * A K5.1 hydrophobic word is five exact amino acids plus one bit derived from the sixth.
 * Raw rolling words retain canonical residues; blacklist lookup chooses raw or projected coordinates.
 * @author Yelan, Keqing
 */
public final class ProteinKmerGeometry{

	public static ProteinKmerGeometry exact(final int k){return make(k, null);}

	public static ProteinKmerGeometry make(final int k, final String bits){
		final String clean=bits==null ? "" : bits.trim();
		if(clean.regionMatches(true, 0, "all:", 0, 4)){return reduced(k, clean.substring(4));}
		if(clean.startsWith("dimer:")){
			final int split=clean.indexOf(':',6);
			if(split<0){throw new IllegalArgumentException("Dimer geometry requires starts and groups");}
			return dimers(k,clean.substring(6,split),clean.substring(split+1));
		}
		if(clean.length()==0 || clean.equalsIgnoreCase("none") || clean.equals("0")){
			return new ProteinKmerGeometry(k, new int[0][], new int[0], "exact"+k);
		}
		final String[] tokens=clean.split(",", -1);
		final int[][] maps=new int[tokens.length][];
		final int[] widths=new int[tokens.length];
		final StringBuilder name=new StringBuilder("exact").append(k);
		for(int i=0; i<tokens.length; i++){
			final String token=tokens[i].trim();
			if(token.length()==0){throw new IllegalArgumentException("Empty amino-acid partition in "+bits);}
			final String partition=partition(token);
			maps[i]=map(partition);
			widths[i]=bits(maps[i]);
			name.append('+').append(partition);
		}
		return new ProteinKmerGeometry(k, maps, widths, name.toString());
	}

	/** Reduces every residue in the word; group order defines the packed symbol codes. */
	public static ProteinKmerGeometry reduced(final int k, final String groups){
		if(groups==null){throw new IllegalArgumentException("Full-word reduction requires an explicit amino-acid partition");}
		final String canonical=partition(groups.trim());
		final int[] map=multiwayMap(canonical, false);
		return new ProteinKmerGeometry(k, map, bits(map), "reduced"+k+"+"+canonical);
	}

	public static ProteinKmerGeometry fromName(final String name){
		if(name!=null && name.startsWith("dimer")){
			final int first=name.indexOf('+'), second=first<0 ? -1 : name.indexOf('+',first+1);
			if(first<=5 || second<0){throw new IllegalArgumentException("Invalid dimer geometry: "+name);}
			return dimers(Integer.parseInt(name.substring(5,first)),name.substring(first+1,second),name.substring(second+1));
		}
		if(name!=null && name.startsWith("reduced")){
			final int split=name.indexOf('+');
			if(split<=7 || split==name.length()-1){throw new IllegalArgumentException("Invalid reduced protein geometry: "+name);}
			return reduced(Integer.parseInt(name.substring(7, split)), name.substring(split+1));
		}
		if(name==null || !name.startsWith("exact")){throw new IllegalArgumentException("Invalid protein k-mer geometry: "+name);}
		final int split=name.indexOf('+');
		final int end=split<0 ? name.length() : split;
		final int k=Integer.parseInt(name.substring(5,end));
		if(split==name.length()-1){throw new IllegalArgumentException("Missing amino-acid partition in "+name);}
		if(split<0){return exact(k);}
		final String suffix=name.substring(split+1).replace('+', ',').trim();
		if(suffix.regionMatches(true, 0, "all:", 0, 4) || suffix.startsWith("dimer:")){
			throw new IllegalArgumentException("Stored projected geometry must use its own identity, not an exact-prefix name: "+name);
		}
		return make(k, suffix);
	}

	private ProteinKmerGeometry(final int exact_, final int[][] maps_, final int[] widths_, final String name_){
		if(exact_<1 || exact_>12 || maps_==null || widths_==null || maps_.length!=widths_.length ||
			exact_+maps_.length>12 || maps_.length>4){
			throw new IllegalArgumentException("Protein k-mer geometry exceeds supported word width");
		}
		int bits=0;
		for(int i=0; i<widths_.length; i++){
			if(maps_[i]==null || widths_[i]<1 || widths_[i]>5){throw new IllegalArgumentException("Invalid reduced amino-acid partition");}
			bits+=widths_[i];
		}
		if(5*exact_+bits>62){throw new IllegalArgumentException("Protein k-mer geometry exceeds supported word width");}
		k=exact_; exact=exact_; maps=maps_; widths=widths_; name=name_; span=exact+maps.length; suffixBits=bits;
		fullMap=null; fullPairs=null; symbolBits=0; dimerPairs=null; dimerShifts=null;
		suffixRawBits=5*maps.length; suffixRawMask=(1<<suffixRawBits)-1; suffixTable=makeSuffixTable(maps, widths, suffixBits);
		keyBits=5*exact+suffixBits;
		exactMask=(1L<<(5*exact))-1;
		rollingMask=(1L<<(5*span))-1;
	}

	private ProteinKmerGeometry(final int k_, final int[] map_, final int bits_, final String name_){
		if(k_<1 || k_>12 || map_==null || map_.length!=32 || bits_<1 || bits_>5 || k_*bits_>62){
			throw new IllegalArgumentException("Reduced protein geometry must fit the raw and projected word widths");
		}
		k=k_; exact=0; span=k_; name=name_; fullMap=map_; symbolBits=bits_; keyBits=k_*bits_;
		dimerPairs=null; dimerShifts=null;
		maps=new int[0][]; widths=new int[0]; suffixBits=0; suffixRawBits=0; suffixRawMask=0;
		suffixTable=null; exactMask=0; rollingMask=(1L<<(5*k_))-1;
		fullPairs=new int[1<<10];
		for(int raw=0; raw<fullPairs.length; raw++){
			fullPairs[raw]=(map_[raw>>>5]<<bits_)|map_[raw&31];
		}
	}

	/** Ordered adjacent-pair classes at explicit one-based starts within a raw K5 or K6 word. */
	public static ProteinKmerGeometry dimers(final int k, final String starts, final String groups){
		if((k!=5 && k!=6) || starts==null || groups==null){
			throw new IllegalArgumentException("Dimer geometry requires K5 or K6, starts and a complete partition");
		}
		final String[] positions=starts.split(",",-1), partitions=groups.split("/",-1);
		if(positions.length<1 || positions.length>=k || partitions.length<2 || partitions.length>400){
			throw new IllegalArgumentException("Invalid dimer start or class count");
		}
		final int[] shifts=new int[positions.length], pairs=new int[1024];
		Arrays.fill(pairs,-1);
		int previous=0, covered=0;
		for(int i=0;i<positions.length;i++){
			final int start=Integer.parseInt(positions[i]);
			if(!positions[i].equals(Integer.toString(start)) || start<=previous || start>=k){
				throw new IllegalArgumentException("Dimer starts must be canonical, increasing positions in1..k-1");
			}
			shifts[i]=5*(k-start-1); covered|=3<<(start-1); previous=start;
		}
		if(covered!=(1<<k)-1){throw new IllegalArgumentException("Dimer starts must cover every residue in the word");}
		int count=0;
		for(int group=0;group<partitions.length;group++){
			for(String pair:partitions[group].split(",",-1)){
				if(pair.length()!=2){throw new IllegalArgumentException("Dimer groups require comma-separated two-residue symbols");}
				final int a=canonicalResidue(pair.charAt(0)), b=canonicalResidue(pair.charAt(1)), raw=(a<<5)|b;
				if(pairs[raw]>=0){throw new IllegalArgumentException("Duplicate dimer: "+pair);}
				pairs[raw]=group; count++;
			}
		}
		if(count!=400){throw new IllegalArgumentException("Dimer partition must cover all400 ordered canonical pairs");}
		final int width=32-Integer.numberOfLeadingZeros(partitions.length-1);
		return new ProteinKmerGeometry(k,pairs,shifts,width,"dimer"+k+"+"+starts+"+"+groups);
	}

	private static int canonicalResidue(final char symbol){
		final int aa=symbol<AminoAcid.acidToNumberNoStops.length ? AminoAcid.acidToNumberNoStops[symbol] : -1;
		if(aa<0 || aa>=20 || AminoAcid.numberToAcid[aa]!=symbol){
			throw new IllegalArgumentException("Dimer partitions require canonical uppercase residues");
		}
		return aa;
	}

	private ProteinKmerGeometry(final int k_, final int[] pairs_, final int[] shifts_, final int bits_, final String name_){
		k=k_; exact=0; span=k_; name=name_; symbolBits=bits_; keyBits=bits_*shifts_.length;
		if(keyBits>60){throw new IllegalArgumentException("Dimer projection must leave space for positional tags");}
		dimerPairs=pairs_; dimerShifts=shifts_; fullMap=null; fullPairs=null;
		maps=new int[0][]; widths=new int[0]; suffixBits=0; suffixRawBits=0; suffixRawMask=0;
		suffixTable=null; exactMask=0; rollingMask=(1L<<(5*k_))-1;
	}

	/** Projects a raw five-bit-per-residue word over {@link #span()} residues into the configured hash word. */
	public long word(final long rolling){
		if(dimerPairs!=null){
			long word=0;
			for(int shift:dimerShifts){
				final int symbol=dimerPairs[(int)((rolling>>>shift)&1023)];
				if(symbol<0){throw new IllegalArgumentException("Dimer word contains a noncanonical raw residue");}
				word=(word<<symbolBits)|symbol;
			}
			assert(word>=0 && word<(1L<<keyBits)) : "Dimer classes must retain their full packed width";
			return word;
		}
		if(fullMap!=null){
			long raw=rolling, word=0;
			int remaining=k, shift=0;
			for(; remaining>=2; remaining-=2, shift+=2*symbolBits){
				word|=((long)fullPairs[(int)(raw&1023)])<<shift; raw>>>=10;
			}
			if(remaining>0){word|=((long)fullMap[(int)(raw&31)])<<shift;}
			assert(word>=0 && word<(1L<<keyBits)) : "Full-word reduction must fit k times its symbol width before hashing";
			return word;
		}
		if(suffixBits==0){return rolling;}
		final long word=((rolling>>>suffixRawBits)<<suffixBits)|suffixTable[(int)(rolling&suffixRawMask)];
		assert(word>=0 && word<(1L<<keyBits)) : "Projected protein word must fit its declared bit width";
		return word;
	}

	/** Returns the exact prefix word before suffix bits are added. */
	public long exactPrefix(final long rolling){
		if(fullMap!=null || dimerPairs!=null){throw new IllegalArgumentException("A fully reduced word has no exact prefix; use the raw pentamer for blacklist checks");}
		return (rolling>>>suffixRawBits)&exactMask;
	}

	/** Returns the first raw amino-acid prefix for exact blacklist checks. */
	public long rawPrefix(final long rolling, final int residues){
		if(residues<1 || residues>span){throw new IllegalArgumentException("Raw prefix must be within the protein-word span");}
		return rolling>>>(5*(span-residues));
	}

	/** Returns the first raw pentamer in this word for exact five-residue blacklist checks. */
	public long pentamerPrefix(final long rolling){
		if(span<5){throw new IllegalArgumentException("Pentamer blacklist requires a five-residue raw span");}
		return rawPrefix(rolling, 5);
	}

	/** Base word length; suffix-only geometries consume their additional residues separately. */
	public int k(){return k;}
	/** Number of unreduced prefix residues; zero for a fully reduced word. */
	public int exact(){return exact;}
	public int span(){return span;}
	public long rollingMask(){return rollingMask;}
	public String name(){return name;}
	public boolean isExact(){return fullMap==null && dimerPairs==null && maps.length==0;}
	public boolean isExact(final int k){return isExact() && exact==k;}
	/** True only when the same explicit amino-acid partition reduces every residue. */
	public boolean isReduced(){return fullMap!=null;}
	public int suffixBits(){return suffixBits;}

	@Override
	public String toString(){return name;}

	@Override
	public int hashCode(){return name.hashCode();}

	@Override
	public boolean equals(Object other){
		if(this==other){return true;}
		if(!(other instanceof ProteinKmerGeometry)){return false;}
		final ProteinKmerGeometry o=(ProteinKmerGeometry)other;
		if(k!=o.k || exact!=o.exact || !name.equals(o.name) || maps.length!=o.maps.length ||
				!Arrays.equals(widths, o.widths) || !Arrays.equals(fullMap, o.fullMap)){return false;}
		for(int i=0; i<maps.length; i++){if(!Arrays.equals(maps[i], o.maps[i])){return false;}}
		return true;
	}

	private static String partition(final String token){
		final String lower=token.toLowerCase(Locale.ROOT);
		if(lower.equals("hydro") || lower.equals("hydro10") || lower.equals("kd10")){return HYDRO_KD_HIGH10;}
		if(lower.equals("blosum") || lower.equals("blosum10") || lower.equals("blosum2")){return BLOSUM62_HIGH10;}
		if(lower.equals("blosum4") || lower.equals("blosum5x4")){return BLOSUM62_5X4;}
		return token.toUpperCase(Locale.ROOT);
	}

	private static int[] map(final String partition){
		return partition.indexOf('/')>=0 ? multiwayMap(partition) : binaryMap(partition);
	}

	private static int[] binaryMap(final String high){
		final int[] map=new int[32];
		final boolean[] seen=new boolean[20];
		for(int i=0; i<high.length(); i++){
			final char symbol=high.charAt(i);
			final int aa=symbol<AminoAcid.acidToNumberNoStops.length ? AminoAcid.acidToNumberNoStops[symbol] : -1;
			if(aa<0 || aa>=20 || AminoAcid.numberToAcid[aa]!=symbol || seen[aa]){
				throw new IllegalArgumentException("Amino-acid partition must list unique canonical uppercase residues: "+high);
			}
			seen[aa]=true; map[aa]=1;
		}
		return map;
	}

	private static int[] multiwayMap(final String partition){
		return multiwayMap(partition, true);
	}

	private static int[] multiwayMap(final String partition, final boolean suffix){
		final String[] groups=partition.split("/", -1);
		if(groups.length<2 || groups.length>(suffix ? 32 : 20) || (suffix && Integer.bitCount(groups.length)!=1)){
			throw new IllegalArgumentException((suffix ? "Suffix partitions need 2, 4, 8, 16 or 32 groups: " :
				"Full-word partitions need 2..20 groups: ")+partition);
		}
		final int[] map=new int[32];
		final boolean[] seen=new boolean[20];
		int count=0;
		for(int group=0; group<groups.length; group++){
			if(groups[group].length()<1){throw new IllegalArgumentException("Empty amino-acid partition group: "+partition);}
			for(int i=0; i<groups[group].length(); i++){
				final char symbol=groups[group].charAt(i);
				final int aa=symbol<AminoAcid.acidToNumberNoStops.length ? AminoAcid.acidToNumberNoStops[symbol] : -1;
				if(aa<0 || aa>=20 || AminoAcid.numberToAcid[aa]!=symbol || seen[aa]){
					throw new IllegalArgumentException("Amino-acid partition must list each canonical uppercase residue once: "+partition);
				}
				seen[aa]=true; map[aa]=group; count++;
			}
		}
		if(count!=20){throw new IllegalArgumentException("Multiway partitions must cover all 20 canonical residues: "+partition);}
		return map;
	}

	private static int bits(final int[] map){
		int max=0;
		for(int i=0; i<20; i++){max=Math.max(max, map[i]);}
		assert(max>=1) : "Reduced amino-acid partitions must use at least two symbols";
		return 32-Integer.numberOfLeadingZeros(max);
	}

	private static int[] makeSuffixTable(final int[][] maps, final int[] widths, final int suffixBits){
		final int[] table=new int[1<<(5*maps.length)];
		for(int raw=0; raw<table.length; raw++){
			int word=0, shift=suffixBits;
			for(int i=0; i<maps.length; i++){
				final int aa=(raw>>>(5*(maps.length-i-1)))&31;
				shift-=widths[i];
				word|=maps[i][aa]<<shift;
			}
			assert(word>=0 && word<(1<<suffixBits)) : "Suffix reduction must fit the declared suffix bits";
			table[raw]=word;
		}
		return table;
	}

	private static final String HYDRO_KD_HIGH10="IVLFCMAGTS";
	private static final String BLOSUM62_HIGH10="ACILMFTWYV";
	private static final String BLOSUM62_5X4="AGPST/RDQEK/NHFWY/CILMV";
	private final int k, exact, span, keyBits, suffixBits, suffixRawBits, suffixRawMask, symbolBits;
	private final long exactMask, rollingMask;
	private final int[][] maps;
	private final int[] widths;
	private final int[] suffixTable;
	private final int[] fullMap, fullPairs;
	private final int[] dimerPairs, dimerShifts;
	private final String name;
}
