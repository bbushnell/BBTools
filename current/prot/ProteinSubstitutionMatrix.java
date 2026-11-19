package prot;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import dna.AminoAcid;
import dna.Data;
import fileIO.ByteFile;
import parse.LineParser1;
import shared.Resources;

/** Immutable integer substitution scores and their linear-gap units.
 * Canonical residues use BBTools codes0..19; ambiguity retains the scaled
 * BLOSUM62 X policy. A portable binding carries the actual scores, not a path
 * or a mutable process-wide selection. Statistical E-value parameters are not
 * inferred for custom matrices.
 * @author Brian, Keqing
 */
public final class ProteinSubstitutionMatrix{

	public static final ProteinSubstitutionMatrix BLOSUM62=standard();
	/** Case-sensitive release name and installed resource selected for new work. */
	public static final String POSSUM55_NAME="POSSUM55";
	public static final String DEFAULT_NAME=POSSUM55_NAME;
	public static final String POSSUM55_RESOURCE="bbgap/POSSUM55.tsv";
	/** Verified name-only successor of the accepted R5 matrix; all scores and units are retained. */
	public static final String POSSUM55_FILE_SHA80="474b6013e399137f0fa1";
	public static final String POSSUM55_BINDING_SHA80="21a59255f7b0dfa77945";
	/** Explicit legacy selection remains available; it is not the construction default. */
	public static final String CUSTOM40_RESOURCE="protein_substitution_custom40.tsv";
	public static final String CUSTOM40_FILE_SHA80="4f1fe89f35bd8ee21f93";
	public static final String CUSTOM40_BINDING_SHA80="f3ddb77a31c0b04e56dc";
	/** Public expansion accepted for the optional native description row; not part of score identity. */
	public static final String POSSUM55_DESCRIPTION="Prokaryotic Ortholog Sequence SUbstitution Matrix";
	/** Detaches validated canonical scores and derives the fixed ambiguity policy and identity. */
	private ProteinSubstitutionMatrix(String name_,int scale_,int gap_,int[][] standard){
		name=name_;scale=scale_;linearGap=gap_;
		scores=new int[22][22];
		for(int[] row:scores){Arrays.fill(row,Integer.MIN_VALUE/4);}
		for(int a=0;a<22;a++){
			if(a==20){continue;}
			for(int b=0;b<22;b++){
				if(b==20){continue;}
				scores[a][b]=a<20 && b<20 ? standard[a][b] : scale*Blosum62.score((byte)a,(byte)b);
			}
		}
		final StringBuilder text=new StringBuilder(name).append(':').append(scale).append(':').append(linearGap);
		for(int a=0;a<20;a++){for(int b=0;b<20;b++){text.append(':').append(scores[a][b]);}}
		binding=text.toString();sha80=DigestSuffix.bytes(binding.getBytes(StandardCharsets.US_ASCII));
		assert(scores.length==22 && binding.length()>400) : "A matrix binding must retain every canonical pair and its score units";
	}

	/** Copies a complete symmetric20x20 table in encoded residue order.
	 * Bounds preserve HbmPositionModel's fixed-point DP overflow contract.
	 */
	public static ProteinSubstitutionMatrix fromScores(String name,int scale,int linearGap,int[][] standard){
		require(name!=null && (POSSUM55_NAME.equals(name) ||
			name.matches("[a-z][a-z0-9_]*") && !name.equalsIgnoreCase(POSSUM55_NAME)),
			"Matrix name must be a lowercase identifier or exactly POSSUM55; lowercase aliases are unsupported");
		require(scale>0 && scale<=16 && linearGap>0 && linearGap<=100,
			"Matrix scale must be1..16 and linear gap1..100 to preserve bounded integer DP arithmetic");
		require(standard!=null && standard.length==20,"Matrix requires exactly20 canonical residue rows");
		for(int a=0;a<20;a++){
			require(standard[a]!=null && standard[a].length==20,"Matrix row "+a+" must contain20 scores");
		}
		for(int a=0;a<20;a++){
			require(standard[a][a]>0,"Matrix diagonal must be positive for self-score normalization: code="+a);
			for(int b=0;b<20;b++){
				final int s=standard[a][b];
				require(s>=-100 && s<=100,"Matrix score outside DP bound[-100,100]: "+a+","+b+"="+s);
				require(s==standard[b][a],"Substitution matrix must be symmetric: "+a+","+b);
				if("blosum62".equals(name)){
					require(scale==1 && linearGap==4 && s==Blosum62.score((byte)a,(byte)b),
						"The blosum62 name is reserved for the unchanged built-in scores and gap4");
				}
			}
		}
		return new ProteinSubstitutionMatrix(name,scale,linearGap,standard);
	}

	/** Names admitted by construction, package loading and optional CLI assertions.
	 * Serialization remains generic; admitting a name never supplies replacement scores.
	 */
	public static boolean supportedSelector(String name){
		return "blosum62".equals(name) || "custom40".equals(name) || "custom60".equals(name)
			|| POSSUM55_NAME.equals(name)
			|| "pairwise40".equals(name) || "blend25".equals(name) || "blend50".equals(name)
			|| "blend75".equals(name) || "consensus25to40".equals(name)
			|| "consensus25to45".equals(name) || "consensus0to35".equals(name);
	}

	/** Resolves the public selector. POSSUM55 and explicitly selected custom40 have
	 * installed accepted resources. Other external matrices need a file; explicit
	 * files retain their own binding and are never silently replaced by the default.
	 */
	public static ProteinSubstitutionMatrix load(String name,String path){
		require(name!=null,"A matrix selector is required");
		if("blosum62".equals(name)){
			require(path==null || path.isEmpty(),"Built-in blosum62 does not accept a replacement matrixfile");
			return BLOSUM62;
		}
		require(supportedSelector(name),"Unknown matrix selector: "+name);
		if(path==null || path.isEmpty()){
			if(POSSUM55_NAME.equals(name)){return installedMatrix(true);}
			if("custom40".equals(name)){return installedMatrix(false);}
		}
		require(path!=null && !path.isEmpty(),"Custom matrix selection requires matrixfile");
		final ProteinSubstitutionMatrix matrix=read(path);
		require(name.equals(matrix.name),"Matrix selector/file name mismatch: requested="+name+" file="+matrix.name);
		return matrix;
	}

	/** New construction chooses POSSUM55; historical readers use their stored binding. */
	public static ProteinSubstitutionMatrix defaultMatrix(){return load(DEFAULT_NAME,null);}
	/** Actual installed default file, for producer input tracking and output-alias protection. */
	public static String defaultResourcePath(){
		return installedResourcePath(POSSUM55_RESOURCE, POSSUM55_NAME);
	}

	/** The historical override applies only to an explicit custom40 selection. */
	private static String custom40ResourcePath(){
		final String override=System.getProperty("bbgap.custom40.resource");
		if(override!=null){
			require(!override.isEmpty(),"bbgap.custom40.resource must not be empty");
			final java.nio.file.Path selected=java.nio.file.Paths.get(override);
			require(selected.isAbsolute() && selected.equals(selected.normalize()) && java.nio.file.Files.isRegularFile(selected),
				"bbgap.custom40.resource requires an absolute normalized regular disk path");
			return selected.toString();
		}
		return installedResourcePath(CUSTOM40_RESOURCE, "custom40");
	}

	/** Resolves an ordinary installed disk resource without a process exit. */
	private static String installedResourcePath(String resourceName, String name){
		// Ordinary file-backed classpath resources need URI decoding before filesystem use.
		final java.net.URL resource=ProteinSubstitutionMatrix.class.getResource("/"+resourceName);
		if(resource!=null && "file".equals(resource.getProtocol())){
			try{return java.nio.file.Paths.get(resource.toURI()).toAbsolutePath().normalize().toString();}
			catch(java.net.URISyntaxException e){throw new IllegalArgumentException("Invalid "+name+" resource URI",e);}
		}
		final String path=Data.findPath("?"+resourceName,false);
		if(path==null){
			throw new Resources.MissingResourceException("ERROR: Required "+name+" matrix is missing: resources/"+
				resourceName+"\nThis small matrix ships with BBTools. Reinstall a complete release from\n"+
				"https://github.com/bbushnell/BBTools/releases/ and keep its resources directory beside the launchers.\n");
		}
		require(!path.startsWith("jar:"),name+" requires a disk resource; jar-only resources are unsupported");
		final java.nio.file.Path resolved=(path.startsWith("file:") ? java.nio.file.Paths.get(java.net.URI.create(path)) :
			java.nio.file.Paths.get(path)).toAbsolutePath().normalize();
		require(java.nio.file.Files.isRegularFile(resolved),name+" requires an installed disk resource: "+resourceName);
		return resolved.toString();
	}

	/** Loads outside class initialization so missing-resource diagnostics are not
	 * wrapped in ExceptionInInitializerError or permanently poison a lazy holder. */
	private static ProteinSubstitutionMatrix installedMatrix(boolean possum){
		ProteinSubstitutionMatrix matrix=possum ? possum55Matrix : custom40Matrix;
		if(matrix!=null){return matrix;}
		synchronized(ProteinSubstitutionMatrix.class){
			matrix=possum ? possum55Matrix : custom40Matrix;
			if(matrix==null){
				final String name=possum ? POSSUM55_NAME : "custom40";
				final String path=possum ? defaultResourcePath() : custom40ResourcePath();
				final String filePin=possum ? POSSUM55_FILE_SHA80 : CUSTOM40_FILE_SHA80;
				final String bindingPin=possum ? POSSUM55_BINDING_SHA80 : CUSTOM40_BINDING_SHA80;
				require(filePin.equals(DigestSuffix.file(path)),name+" resource byte identity differs");
				matrix=read(path);
				require(name.equals(matrix.name) && matrix.scale==4 && matrix.linearGap==16 &&
					bindingPin.equals(matrix.sha80),name+" scores or units differ from accepted binding");
				if(possum){possum55Matrix=matrix;}else{custom40Matrix=matrix;}
			}
		}
		return matrix;
	}
	private static volatile ProteinSubstitutionMatrix possum55Matrix, custom40Matrix;

	/** Detects native TSV or NCBI text by content. Plain NCBI input needs explicit
	 * units through ProteinSubstitutionMatrixIO; our NCBI exports carry those units.
	 * Imported ambiguity/stop cells never change the established runtime policy.
	 */
	public static ProteinSubstitutionMatrix read(String path){
		return ProteinSubstitutionMatrixIO.read(path, null, 0, 0);
	}

	/** Reads the small version1 TSV: name, scale, gap and400 score/letter/letter/integer
	 * rows. POSSUM55 may include its exact public expansion in one description row.
	 * The description does not change score identity. Provenance tables remain separate;
	 * unknown/duplicate rows fail, and historical files need no new metadata.
	 */
	static ProteinSubstitutionMatrix readNative(String path, ByteFile input, byte[] first){
		assert(input!=null) : "The format dispatcher owns the native input and its close operation";
		final LineParser1 lp=new LineParser1('\t');
		final int[][] table=new int[20][20];final boolean[][] seen=new boolean[20][20];
		String name=null;int scale=0,gap=0,count=0;
		boolean described=false;
		{
			require(Arrays.equals(first,"#bbgap_substitution_matrix\t1".getBytes(StandardCharsets.US_ASCII)),
				"Unexpected substitution matrix schema: "+path);
			for(byte[] line=input.nextLine();line!=null;line=input.nextLine()){
				lp.set(line);
				if(lp.terms()==4 && lp.termEquals("score",0)){
					final int a=residue(lp.parseString(1)),b=residue(lp.parseString(2));
					require(!seen[a][b],"Duplicate matrix score: "+a+","+b);
					table[a][b]=Integer.parseInt(lp.parseString(3));seen[a][b]=true;count++;
				}else{
					require(lp.terms()==2,"Matrix metadata requires exactly two fields");
					if(lp.termEquals("name",0)){require(name==null,"Repeated matrix name");name=lp.parseString(1);}
					else if(lp.termEquals("description",0)){
						require(!described,"Repeated matrix description");
						require(lp.termEquals(POSSUM55_DESCRIPTION,1),"Unsupported matrix description");
						described=true;
					}else if(lp.termEquals("scale",0)){require(scale==0,"Repeated matrix scale");scale=Integer.parseInt(lp.parseString(1));require(scale>0,"Matrix scale must be positive");}
					else if(lp.termEquals("gap",0)){require(gap==0,"Repeated matrix gap");gap=Integer.parseInt(lp.parseString(1));require(gap>0,"Matrix gap must be positive");}
					else{throw new IllegalArgumentException("Unknown matrix metadata: "+lp.parseString(0));}
				}
			}
		}
		require(count==400 && name!=null && scale>0 && gap>0,"Incomplete matrix resource: "+path);
		require(!described || POSSUM55_NAME.equals(name),"Matrix description is reserved for POSSUM55");
		return fromScores(name,scale,gap,table);
	}

	/** Parses the self-contained package-header binding once at load time.
	 * This does not reinterpret a name using a different external resource.
	 */
	public static ProteinSubstitutionMatrix parseBinding(String text){
		require(text!=null,"Missing matrix binding");
		if("blosum62".equals(text)){return BLOSUM62;}
		final String[] parts=text.split(":",-1);
		require(parts.length==403,"Matrix binding requires name, scale, gap and400 scores");
		final int[][] table=new int[20][20];int at=3;
		for(int a=0;a<20;a++){for(int b=0;b<20;b++){table[a][b]=Integer.parseInt(parts[at++]);}}
		return fromScores(parts[0],Integer.parseInt(parts[1]),Integer.parseInt(parts[2]),table);
	}

	/** Compact portable record; caller must retain it with calibrated thresholds. */
	public String binding(){return binding;}
	/** Equality includes all scores and gap units, never a resource's filename alone. */
	public boolean same(ProteinSubstitutionMatrix other){return other!=null && binding.equals(other.binding);}
	/** Canonical or X encoded substitution score; encoded stop20 is invalid. */
	public int score(byte a,byte b){
		assert(a>=0 && a<22 && a!=20 && b>=0 && b<22 && b!=20) : "Matrix scoring requires canonical/X encoded residues";
		return scores[a][b];
	}
	/** Sum of diagonal scores in this matrix's units, with explicit overflow check. */
	public int selfScore(byte[] sequence){
		require(sequence!=null,"Self-score requires an encoded sequence");
		long total=0;
		for(byte a:sequence){require(a>=0 && a<22 && a!=20,"Invalid self-score residue: "+a);total+=scores[a][a];}
		require(total>=Integer.MIN_VALUE && total<=Integer.MAX_VALUE,"Matrix self-score exceeds integer range");
		return (int)total;
	}
	/** Maps an explicit resource letter; ambiguity and encoded stop are never matrix rows. */
	private static int residue(String text){
		require(text.length()==1 && text.charAt(0)>='A' && text.charAt(0)<='Z',"Score residue must be one uppercase canonical letter");
		final int a=AminoAcid.acidToNumberExtended[text.charAt(0)];
		require(a>=0 && a<20,"Matrix resource cannot include ambiguity or stop rows: "+text);return a;
	}
	/** Captures the published built-in values without altering Blosum62 or its callers. */
	private static ProteinSubstitutionMatrix standard(){
		final int[][] table=new int[20][20];
		for(int a=0;a<20;a++){for(int b=0;b<20;b++){table[a][b]=Blosum62.score((byte)a,(byte)b);}}
		return fromScores("blosum62",1,4,table);
	}
	private static void require(boolean valid,String why){if(!valid){throw new IllegalArgumentException(why);}}
	public final String name,sha80;
	public final int scale,linearGap;
	private final int[][] scores;
	private final String binding;
}
