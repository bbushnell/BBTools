package prot;

/** Small, explicit gap-inclusive identity intervals for overlapping recounts. */
public final class IdentityBands{

	private IdentityBands(String[] labels_, double[] min_, double[] max_){
		labels=labels_; min=min_; max=max_;
	}

	public static IdentityBands legacy(){
		return new IdentityBands(new String[]{"<40", "40-60", ">=60"}, new double[]{0, 40, 60}, new double[]{40, 60, 101});
	}

	public static IdentityBands parse(String s){
		if(s==null || s.equalsIgnoreCase("legacy")){return legacy();}
		final String[] parts=s.split(",", -1);
		require(parts.length>0 && parts.length<=30, "At least 1 and at most 30 identity bands are required");
		final String[] labels=new String[parts.length];
		final double[] min=new double[parts.length], max=new double[parts.length];
		for(int i=0; i<parts.length; i++){
			final String part=parts[i];
			final int mid=part.indexOf("to");
			require(mid>0 && mid==part.lastIndexOf("to") && mid+2<part.length(),
				"Custom bands use labels like 25to40: "+part);
			min[i]=Double.parseDouble(part.substring(0, mid));
			max[i]=Double.parseDouble(part.substring(mid+2));
			require(min[i]>=0 && max[i]<=100 && min[i]<max[i],
				"Identity band must satisfy 0<=min<max<=100: "+part);
			labels[i]=part;
			for(int j=0; j<i; j++){require(!labels[i].equals(labels[j]), "Duplicate identity band: "+part);}
		}
		return new IdentityBands(labels, min, max);
	}

	public int mask(double identity){
		require(Double.isFinite(identity) && identity>=0 && identity<=100, "Invalid alignment identity: "+identity);
		int mask=0;
		for(int i=0; i<labels.length; i++){
			if(identity>=min[i] && identity<max[i]){mask|=1<<i;}
		}
		return mask;
	}

	public int exclusiveBand(double identity){
		final int mask=mask(identity);
		require(mask!=0 && (mask&(mask-1))==0, "Expected one identity band, found mask "+mask+" for "+identity);
		return Integer.numberOfTrailingZeros(mask);
	}

	public String label(int band){return labels[band];}
	public int size(){return labels.length;}
	public String recipe(){return String.join(",", labels);}
	public String[] labels(){return labels.clone();}

	private static void require(boolean condition, String message){
		if(!condition){throw new IllegalArgumentException(message);}
	}

	private final String[] labels;
	private final double[] min, max;
}
