package prot;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.HashSet;

import fileIO.ByteFile;
import parse.LineParser1;
import structures.ByteBuilder;
import structures.IntHashMap;

/** Immutable five-scalar member means bound to one loaded family roster.
 * Reads the native FamilyMemberScalars all-mode table without recomputing
 * biological values. Empty physical rows remain undefined and cannot supply
 * values for live families. No consensus or zero-value fallback is allowed.
 * @author Brian Bushnell, Keqing
 */
public final class FamilyMemberScalarMeans{

	public FamilyMemberScalarMeans(FamilyTailCascade.LoadedClassifier classifier,Path table,String fingerprint){
		require(classifier!=null && table!=null && table.isAbsolute(),
			"Member means require a loaded classifier and an absolute table path");
		bound=classifier.config();path=table.normalize();pin=fingerprint;
		final int count=bound.families.length;
		require(count>0 && count<=Integer.MAX_VALUE/COUNT,"Five member means must fit each live family ordinal");
		means=new double[count*COUNT];
		verify();read();verify();
	}

	/** Copy raw means into caller-owned scratch before the existing combined25 normalization. */
	public void copyMeans(int ordinal,double[] out){
		require(ordinal>=0 && ordinal<bound.families.length && out!=null && out.length==COUNT,
			"Combined25 requires five raw member means for a valid live family ordinal");
		System.arraycopy(means,ordinal*COUNT,out,0,COUNT);
	}

	/** Gates may change in a derived view; the index and retained-member roster must not. */
	boolean boundTo(FamilyTailCascade.SearchConfig candidate){
		if(candidate==null || candidate.index!=bound.index || candidate.families.length!=bound.families.length){return false;}
		for(int i=0;i<bound.families.length;i++){
			final FamilyReference a=bound.families[i],b=candidate.families[i];
			if(a.rank!=b.rank || a.members!=b.members || !a.id.equals(b.id)){return false;}
		}
		return true;
	}

	/** Reject changed contents before and after an evaluation using this immutable binding. */
	public void verify(){
		require(pin!=null && pin.matches("[0-9a-f]{20}") && Files.isRegularFile(path) && !Files.isSymbolicLink(path),
			"Member scalar means require a regular file and a lowercase sha80: "+path);
		require(pin.equals(DigestSuffix.file(path.toString())),"Member scalar mean table changed: "+path);
	}
	public Path path(){return path;}
	public String sha80(){return pin;}

	/** Native all-mode output has dense physical ranks; bind only nonempty rows to live ordinals. */
	private void read(){
		final IntHashMap ordinals=new IntHashMap(bound.families.length*2+1);
		for(int i=0;i<bound.families.length;i++){
			final int rank=bound.families[i].rank;
			require(rank>=0 && ordinals.get(rank)<0,"Live family ranks must be unique for member scalar binding");
			ordinals.put(rank,i);
		}
		final HashSet<String> ids=new HashSet<String>();
		final boolean[] seen=new boolean[bound.families.length];
		final ByteFile input=ByteFile.makeByteFile1(path.toString(),false);
		final LineParser1 lp=new LineParser1('\t');
		int rows=0,live=0;
		try{
			require(Arrays.equals(input.nextLine(),header()),"Member means require the native23-column scalar-moment schema");
			for(byte[] line=input.nextLine();line!=null;line=input.nextLine()){
				lp.set(line);
				require(lp.terms()==23,"Member scalar rows require rank/ID/N plus five mean/SD/min/max groups");
				final int rank=lp.parseInt(0);final String id=lp.parseString(1);final long n=lp.parseLong(2);
				require(rank==rows && n>=0 && !id.isEmpty() && ids.add(id),
					"Native member scalar rows must have dense physical ranks, unique identities and nonnegative counts");
				rows++;
				final int ordinal=ordinals.get(rank);
				if(n==0){
					require(ordinal<0,"A live family cannot receive undefined zero-member scalar means: "+id);
					for(int c=3;c<23;c++){require(lp.termEquals("NA",c),"Empty family statistics must remain explicitly undefined: "+id);}
					continue;
				}
				require(ordinal>=0 && !seen[ordinal] && bound.families[ordinal].id.equals(id) && bound.families[ordinal].members==n,
					"Member scalar rank, identity and count must match the exact live classifier roster: "+id);
				for(int scalar=0;scalar<COUNT;scalar++){
					final int c=3+scalar*4;
					final double mean=lp.parseDouble(c),sd=lp.parseDouble(c+1),min=lp.parseDouble(c+2),max=lp.parseDouble(c+3);
					require(Double.isFinite(mean) && Double.isFinite(sd) && Double.isFinite(min) && Double.isFinite(max) &&
						min>=LOW[scalar] && max<=HIGH[scalar] && min<=mean && mean<=max && sd>=0 && sd<=HIGH[scalar]-LOW[scalar],
						"Member scalar moments must be finite and obey native bounds/min<=mean<=max: "+id+" "+NAMES[scalar]);
					require(n!=1 || sd==0 && min==mean && mean==max,"Singleton scalar moments must describe its one member: "+id);
					means[ordinal*COUNT+scalar]=mean;
				}
				seen[ordinal]=true;live++;
			}
		}finally{require(!input.close(),"Member scalar mean table read failed: "+path);}
		require(live==bound.families.length,"Member scalar means must cover every live family; no consensus fallback");
	}

	private static byte[] header(){
		assert(NAMES.length==COUNT && METRICS.length==4) : "FamilyMemberScalars all-mode exports four moments for each native scalar";
		final ByteBuilder bb=new ByteBuilder("family_rank\tfamily_id\tN");
		for(String name:NAMES){for(String metric:METRICS){bb.tab().append(name).append('_').append(metric);}}
		return bb.toBytes();
	}
	private static void require(boolean valid,String why){if(!valid){throw new IllegalArgumentException(why);}}
	private final FamilyTailCascade.SearchConfig bound;
	private final Path path;
	private final String pin;
	private final double[] means;
	private static final int COUNT=ConsensusScalars.SCALAR_COUNT;
	private static final String[] NAMES={"hydrophobicity","net_charge","aromaticity","fraction_small","pI_tail"};
	private static final String[] METRICS={"mean","sd","min","max"};
	/** ConsensusScalars KD values, normalized charge/fractions, and the native pI bracket. */
	private static final double[] LOW={-4.5,-1,0,0,0},HIGH={4.5,1,1,1,20};
}
