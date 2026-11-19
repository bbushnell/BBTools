package prot;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import ddl.DDLFamilySearch.StructuralInputs;
import ddl.ProteinPositionMode;
import fileIO.ByteFile;
import ml.CellNet;
import ml.CellNetParser;
import parse.LineParser1;
import prot.FamilyTailCascade.SearchConfig;
import structures.ByteBuilder;

/** Immutable frozen-model and family-length binding for opt-in classifier ranking.
 * Model state and same-query prediction caches belong to individual workers.
 * @author Keqing
 */
public final class FamilyLengthScore{
	public enum Mode{OFF,MULTIPLY,ADD}

	FamilyLengthScore(SearchConfig config,Path profiles_,String profilePin_,Path model_,String modelPin_,Path lengths_,String lengthPin_) throws Exception{
		this(config,profiles_,profilePin_,model_,modelPin_,lengths_,lengthPin_,
			FamilyCalibration.VerificationMode.SHA80,null);
	}

	FamilyLengthScore(SearchConfig config,Path profiles_,String profilePin_,Path model_,String modelPin_,
			Path lengths_,String lengthPin_,FamilyCalibration.VerificationMode verificationMode_,
			StructuralInputs structuralInputs_) throws Exception{
		require(config!=null && config.top==50 && config.geometry.isExact(5) && config.positional && config.width==1024 &&
			config.features==FamilyTailCascade.RankingFeatures.COMBINED25 && config.cosineMode==FamilyTailCascade.CosineMode.RANKING &&
			config.truth==null && !config.union && !config.candidateOnly && config.ddlReject>=1 && config.ddlAccept>1024 && !config.index.usesPruning(),
			"Released learned classifier ranking requires exact positionalK5, combined25 ranking, strict top50, no index pruning or DDL shortcut");
		require(verificationMode_!=null &&
			(verificationMode_==FamilyCalibration.VerificationMode.STRUCTURAL)==(structuralInputs_!=null),
			"Length ranking requires a structural tracker if and only if verification is structural");
		verificationMode=verificationMode_;structuralInputs=structuralInputs_;
		families=config.families; index=config.index;
		profiles=profiles_.toAbsolutePath().normalize(); profilePin=profilePin_;
		model=model_.toAbsolutePath().normalize(); modelPin=modelPin_;
		lengths=lengths_.toAbsolutePath().normalize(); lengthPin=lengthPin_;
		verify(); require(Files.size(model)<=1000000,"The released rank model is a small frozen network, not an unbounded model load");
		modelLines=ByteFile.toLines(model.toString()); validateModel(CellNetParser.loadInferenceFromLines(modelLines));
		medians=readMedians();
		final char[][] rows=FamilyTailCascade.loadReplacementRows(families,profiles,config.geometry,
			ProteinPositionMode.FIXED,config.partitions,config.overlap,config.blacklist,structuralInputs);
		occupied=new int[rows.length];
		for(int f=0;f<rows.length;f++){for(char value:rows[f]){if(value!=0){occupied[f]++;}} rows[f]=null;}
		verify();
	}

	/** Reject changed stored bytes, including compressed table bytes, before publishing a result. */
	public void verify() throws Exception{
		verifyInput(profiles,profilePin,"Length ranking profiles");
		verifyInput(model,modelPin,"Length ranking model");
		verifyInput(lengths,lengthPin,"Length ranking medians");
		if(structuralInputs!=null){structuralInputs.verify();}
	}
	public ByteBuilder recipe(){
		return new ByteBuilder("length_model\t").append(model.toString()).append("\nlength_model_sha80\t").append(modelPin)
			.append("\nfamily_medians\t").append(lengths.toString()).append("\nfamily_medians_sha80\t").append(lengthPin)
			.append("\nlength_profiles_sha80\t").append(profilePin)
			.append("\nlength_verification\t").append(verificationMode==FamilyCalibration.VerificationMode.STRUCTURAL ?
				"structural\nlength_identity_semantics\tinherited_unverified" :
				"sha80\nlength_identity_semantics\tmeasured_sha80")
			.append("\nlength_features\tln_query_length_over_member_median,matched_family_register_fraction,combined25_cosine\n")
			.append("length_normalization\tmodel_owned_once\nlength_add_weight\t1\nlength_query_range\t100..199\n");
	}
	boolean boundTo(SearchConfig config){return config.index==index && config.families==families;}
	Worker newWorker(){return new Worker();}

	final class Worker{
		Worker(){net=CellNetParser.loadInferenceFromLines(modelLines);validateModel(net);}
		/** A generation denotes one raw query, shared by its three ranking variants. */
		float predict(int length,int ordinal,int hits,double cosine,long generation){
			require(generation>0 && ordinal>=0 && ordinal<families.length,"Prediction cache needs a current query generation and live family ordinal");
			if(stamps[ordinal]!=generation){
				require(length>=100 && length<=199 && hits>0 && occupied[ordinal]>=hits && Double.isFinite(cosine) && cosine>=0 && cosine<=1,
					"Length ranking needs100..199-aa queries and finite features from the bound family profile");
				input[0]=(float)Math.log(length/medians[ordinal]);input[1]=(float)(hits/(double)occupied[ordinal]);input[2]=(float)cosine;
				net.applyInput(input);final float f=net.feedForward();
				require(Float.isFinite(f) && f>=0 && f<=1,"Frozen length model produced a nonfinite or non-sigmoid value");
				predictions[ordinal]=f;stamps[ordinal]=generation;
			}
			return predictions[ordinal];
		}
		final CellNet net;final float[] input=new float[3],predictions=new float[families.length];
		final long[] stamps=new long[families.length];
	}

	private double[] readMedians(){
		final HashMap<String,Integer> ordinals=new HashMap<String,Integer>();
		for(int f=0;f<families.length;f++){require(ordinals.put(families[f].id,f)==null,"Live family IDs must be unique");}
		final double[] out=new double[families.length];final boolean[] seen=new boolean[out.length];
		final ByteFile input=ByteFile.makeByteFile1(lengths.toString(),false);final LineParser1 lp=new LineParser1('\t');int count=0;
		try{
			require(Arrays.equals(input.nextLine(),"family\tmembers\tlength_p10\tlength_median\tlength_p90\tcount_under_200".getBytes(StandardCharsets.US_ASCII)),"Family median table has an unexpected schema");
			for(byte[] line=input.nextLine();line!=null;line=input.nextLine()){
				lp.set(line);require(lp.terms()==6,"Family median rows need all six census fields");
				final Integer ordinal=ordinals.get(lp.parseString(0));
				require(ordinal!=null && !seen[ordinal] && lp.parseLong(1)==families[ordinal].members,"Family medians must match live IDs and member counts exactly once");
				final double median=lp.parseDouble(3);
				require(Double.isFinite(median) && median>=100 && median<=200,"A live family needs a finite retained-member median in100..200; no imputation");
				out[ordinal]=median;seen[ordinal]=true;count++;
			}
		}finally{require(!input.close(),"Family median table read failed");}
		require(count==families.length,"Family median table is incomplete");return out;
	}
	private static void validateModel(CellNet net){
		require(net.numInputs()==3 && net.numOutputs()==1 && net.dims.length<=8,"Frozen rank network needs three inputs, one output and bounded topology");
		for(int width:net.dims){require(width>0 && width<=32,"Frozen rank network exceeds the released small topology");}
	}
	private static void pin(Path path,String expected) throws Exception{
		require(expected!=null && expected.matches("[0-9a-f]{20}") && Files.isRegularFile(path) && !Files.isSymbolicLink(path) && expected.equals(DigestSuffix.file(path.toString())),"Length-ranking binding fingerprint differs: "+path);
	}
	private void verifyInput(Path path,String expected,String label)throws Exception{
		if(structuralInputs==null){pin(path,expected);}
		else{
			require(expected!=null && (expected.matches("[0-9a-f]{20}") || expected.equals("unmeasured")),
				label+" needs an inherited lowercase sha80 or unmeasured label");
			structuralInputs.track(path);
		}
	}
	private static void require(boolean valid,String why){if(!valid){throw new IllegalArgumentException(why);}}
	private final FamilyReference[] families;
	private final ddl.DDLFamilySearch.Index index;
	private final FamilyCalibration.VerificationMode verificationMode;
	private final StructuralInputs structuralInputs;
	private final Path profiles,model,lengths;
	private final String profilePin,modelPin,lengthPin;
	private final double[] medians;
	private final int[] occupied;
	private final ArrayList<byte[]> modelLines;
}
