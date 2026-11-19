package prot;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;

import ddl.DDLFamilySearch.StructuralInputs;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import structures.ByteBuilder;

/** Pins a saved floor vector to the actual parent classifier recipe, absolute
 * margin and fitting-cohort manifest. This context loads its own parent, so a
 * caller cannot supply an unrelated resident classifier with a matching label.
 * The fitting driver must independently verify cohort members and replay results;
 * this codec checks provenance and counts, not observed scientific error rates.
 * @author Brian, Keqing
 */
public final class FamilyHitFloorBinding{

	public FamilyHitFloorBinding(final Path parent,final String parentPin,final double margin_,
			final Path cohorts_,final String cohortPin_,final int threads)throws Exception{
		this(parent,parentPin,margin_,cohorts_,cohortPin_,threads,FamilyCalibration.VerificationMode.SHA80);
	}

	public FamilyHitFloorBinding(final Path parent,final String parentPin,final double margin_,
			final Path cohorts_,final String cohortPin_,final int threads,
			final FamilyCalibration.VerificationMode verificationMode_)throws Exception{
		require(parent!=null && cohorts_!=null && threads>0
				&& Double.isFinite(margin_) && margin_>=0 && verificationMode_!=null,
			"A floor binding needs declared parent/cohort files, positive workers and a finite absolute margin");
		structuralInputs=verificationMode_==FamilyCalibration.VerificationMode.STRUCTURAL ? new StructuralInputs() : null;
		parentLabel=verifyInput(parent,parentPin,"Parent recipe");
		cohorts=cohorts_.toAbsolutePath().normalize(); margin=margin_;
		cohortLabel=verifyInput(cohorts,cohortPin_,"Fitting cohort manifest");
		pack=new FamilyScoringPackage(parent,verificationMode_);
		require(structuralInputs!=null || pack.fingerprint.equals(parentLabel),
			"Parent recipe changed while opening its floor binding");
		base=pack.load(threads).withEvaluationMarginScale(margin);
		// Validate the released short-family geometry before reading or writing any vector.
		base.withFamilyHitFloors(base.newFamilyHitFloorProbe());
		verify();
	}

	/** One shared loaded parent for fitting and all exact replays; no reload per vector. */
	public FamilyTailCascade.LoadedClassifier parent(){return base;}
	public FamilyHitFloorVector probe(){return base.newFamilyHitFloorProbe();}
	/** Native fitting reads the exact cohort manifest already bound to this context. */
	Path fittingManifest(){return cohorts;}

	/** Write a requested final raw TSV. The execution controller owns size-appropriate
	 * verified compression. No intermediate matrix or sequence file is emitted.
	 * Initial winner bins are the actual floor-one positive measurement, not a later
	 * replay. Unsupported families must retain six in the selected vector. */
	public void write(final Path output,final FamilyHitFloorVector vector,final long positiveQueries,
			final long[] initialWinners)throws Exception{
		require(output!=null && output.toString().endsWith(".tsv")
				&& !Files.exists(output) && !Files.isSymbolicLink(output),
			"Write a new final raw .tsv floor table; existing artifacts and compression are not overwritten here");
		require(vector!=null && vector.hasInitialSupport() && vector.boundToPolicy(base.config()),
			"Only a supported vector bound to this actual parent policy can be serialized");
		require(!parentLabel.equals("unmeasured") && !cohortLabel.equals("unmeasured"),
			"Structural floor-table writes require explicit declared parent and cohort labels");
		final FamilyHitFloorVector initial=probe().withInitialSupport(positiveQueries,initialWinners);
		require(Arrays.equals(initial.copyInitialSupport(),vector.copyInitialSupport()),
			"Serialized initial winner bins must reproduce the vector's frozen fitting support");
		base.withFamilyHitFloors(vector);
		verify();
		final FamilyReference[] families=base.config().families;
		final ByteBuilder line=new ByteBuilder();
		final ByteStreamWriter writer=new ByteStreamWriter(output.toString(),false,false,false);
		writer.start();
		try{
			line.append(MAGIC).nl().append("parent_sha80\t").append(parentLabel).nl()
				.append("cohorts_sha80\t").append(cohortLabel).nl()
				.append("margin_scale\t").append(Double.toString(margin)).nl()
				.append("positive_queries\t").append(positiveQueries).nl()
				.append(HEADER).nl();
			writer.print(line);
			for(int f=0; f<families.length; f++){
				final FamilyReference family=families[f];
				line.clear().append(family.rank).tab().append(family.id).tab().append(family.members).tab()
					.append(family.residues).tab().append(vector.floor(f));
				for(int hit=0; hit<BINS; hit++){line.tab().append(initialWinners[f*BINS+hit]);}
				writer.print(line.nl());
			}
		}finally{require(!writer.poisonAndWait(),
			"Final floor-table write failed; the partial output is not a completed vector");}
		verify();
	}

	/** Load one complete pinned table and attach it to this context's actual parent.
	 * Physical rank gaps are retained; rows must follow the exact live reference order. */
	public FamilyTailCascade.LoadedClassifier load(final Path table,final String tablePin)throws Exception{
		verify(); verifyInput(table,tablePin,"Floor table");
		final FamilyReference[] families=base.config().families;
		final byte[] floors=new byte[families.length];
		final long[] winners=new long[Math.multiplyExact(families.length,BINS)];
		final ByteFile input=ByteFile.makeByteFile1(table.toString(),false);
		final LineParser1 lp=new LineParser1('\t');
		long queries=0; int rows=0;
		try{
			require(Arrays.equals(input.nextLine(),ascii(MAGIC)),"Unsupported family-floor table schema");
			metadata(input,lp,"parent_sha80",parentLabel);
			metadata(input,lp,"cohorts_sha80",cohortLabel);
			metadata(input,lp,"margin_scale",Double.toString(margin));
			final byte[] denominator=input.nextLine();
			require(denominator!=null,"Missing initial fitting-positive denominator");
			lp.set(denominator);
			require(lp.terms()==2 && lp.termEquals("positive_queries",0),
				"Expected the initial fitting-positive query denominator");
			queries=count(lp,denominator,1,Long.MAX_VALUE);
			require(queries>0,"Initial fitting positives require a positive observed denominator");
			require(Arrays.equals(input.nextLine(),ascii(HEADER)),
				"Floor rows require exact rank, identity, census, floor and seven initial winner bins");
			for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
				require(rows<families.length,"Floor table contains extra families outside the loaded live roster");
				lp.set(line); final FamilyReference family=families[rows];
				require(lp.terms()==5+BINS && count(lp,line,0,Integer.MAX_VALUE)==family.rank
					&& lp.termEquals(family.id,1) &&
					count(lp,line,2,Long.MAX_VALUE)==family.members && count(lp,line,3,Long.MAX_VALUE)==family.residues,
					"Floor rows must match every live reference's physical rank, identity and member census in order");
				final int floor=(int)count(lp,line,4,6);
				require(floor>=1 && floor<=6,"Serialized family floors must remain within the released one-through-six range");
				floors[rows]=(byte)floor;
				for(int hit=0; hit<BINS; hit++){winners[rows*BINS+hit]=count(lp,line,5+hit,Long.MAX_VALUE);}
				rows++;
			}
		}finally{require(!input.close(),"Floor-table read failed; no partial roster may be attached");}
		require(rows==families.length,"Floor table omits live families; no fixed-floor imputation is allowed");
		final FamilyHitFloorVector vector=probe().withInitialSupport(queries,winners).withFloors(floors);
		verifyInput(table,tablePin,"Floor table"); verify();
		return base.withFamilyHitFloors(vector);
	}

	/** The package rechecks all bound classifier components; the cohort manifest
	 * fingerprint identifies fitting inputs whose individual pins the driver checks. */
	public void verify()throws Exception{
		pack.verify();verifyInput(cohorts,cohortLabel,"Fitting cohort manifest");
		if(structuralInputs!=null){structuralInputs.verify();}
	}
	/** Parse unsigned decimal counts with checked accumulation. The ordinary native
	 * integer parser does not detect overflow, which could turn corrupt counts into
	 * plausible positive measurements before the conservation checks see them. */
	private static long count(final LineParser1 lp,final byte[] line,final int term,final long max){
		assert(max>=0) : "Floor-table count bounds must be nonnegative to prevent overflow";
		final int length=lp.setBounds(term);
		if(length==0){throw new IllegalArgumentException("Empty floor-table count in field "+term);}
		long value=0;
		for(int i=lp.a(); i<lp.b(); i++){
			final int digit=line[i]-'0';
			if(digit<0 || digit>9 || value>max/10 || (value==max/10 && digit>max%10)){
				throw new IllegalArgumentException("Floor-table field "+term+" must be unsigned decimal within 0.."+max);
			}
			value=value*10+digit;
		}
		return value;
	}
	private static void metadata(final ByteFile input,final LineParser1 lp,final String key,final String expected){
		final byte[] line=input.nextLine(); require(line!=null,"Incomplete floor-table provenance: "+key);
		lp.set(line);
		require(lp.terms()==2 && lp.termEquals(key,0) && lp.termEquals(expected,1),
			"Floor-table provenance differs from the bound parent, margin or fitting cohort: "+key);
	}
	private static byte[] ascii(final String text){return text.getBytes(StandardCharsets.US_ASCII);}
	private String verifyInput(final Path path,final String expected,final String label)throws Exception{
		require(path!=null && expected!=null &&
			(expected.matches("[0-9a-f]{20}") || structuralInputs!=null && expected.equals("unmeasured")),
			label+" requires a regular file and declared lowercase sha80 or structural unmeasured label");
		if(structuralInputs!=null){
			structuralInputs.track(path);
		}else{
			require(Files.isRegularFile(path) && !Files.isSymbolicLink(path),
				label+" requires a regular non-symlink file");
			require(expected.equals(DigestSuffix.file(path.toString())),label+" changed: "+path);
		}
		return expected;
	}
	private static void require(final boolean valid,final String why){if(!valid){throw new IllegalArgumentException(why);}}
	private final FamilyScoringPackage pack;
	private final FamilyTailCascade.LoadedClassifier base;
	private final StructuralInputs structuralInputs;
	private final Path cohorts;
	private final String parentLabel,cohortLabel;
	private final double margin;
	private static final int BINS=FamilyHitFloorSelector.BINS;
	private static final String MAGIC="#family_hit_floors\t1";
	private static final String HEADER="family_rank\tfamily_id\tmembers\tresidues\tfloor"
		+ "\tinitial_hits0\tinitial_hits1\tinitial_hits2\tinitial_hits3"
		+ "\tinitial_hits4\tinitial_hits5\tinitial_hits6plus";
}
