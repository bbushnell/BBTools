package fun;

public class ProbShared {

	public static void main(String args[]){
		int k=Integer.parseInt(args[0]);
		int len1=Integer.parseInt(args[1]);
		int len2=Integer.parseInt(args[2]);

		System.out.println("Cardinality 1: "+cardinality(k, len1));
		System.out.println("Cardinality 2: "+cardinality(k, len2));
		System.out.println("Probability:   "+probIntersect(k, len1, len2));
		
	}
	
	static int cardinality(int k, int seqLength){
		double space=Math.pow(4, k);
		int kmers=seqLength-k+1;
		double unique=0;
		for(int i=0; i<kmers; i++){
			double prob=(space-unique)/space;
			unique+=prob;
		}
		return (int)Math.round(unique);
	}

	static double probIntersect(int k, int len1, int len2){
		int card1=cardinality(k, len1);
		int card2=cardinality(k, len2);
		double space=Math.pow(4, k);
		double cumulativeProbUnshared=1;
		for(int i=0; i<card1; i++){
			double probShared=card2/space;
			double probUnshared=1-probShared;
			space-=probUnshared;
			cumulativeProbUnshared*=probUnshared;
		}
		return 1-cumulativeProbUnshared;
	}
	
}
