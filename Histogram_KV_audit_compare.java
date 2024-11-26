package kv_whole;//package kv_whole;

import com.google.common.collect.Ordering;
import com.google.common.collect.TreeMultimap;
import keyvalue.data.Processor;
import org.apache.commons.math3.distribution.BetaDistribution;

import java.text.DecimalFormat;
import java.util.*;
import java.util.stream.IntStream;
//import javafx.util.Pair;

//todo: hier 2 bit is VKV-Auditor ue formatted is HKV-Auditor output

public class Histogram_KV_audit_compare {
	int iter = 10000000;
	// 0: AppData    1: JData     2: gaussian    3: powerlaw     4: amazon    5: cloting
	int dataset_index = 4;   // 0-5
	int[] datasets =  {11,    12,    21,    22,    31,      32};
	int dataset = datasets[dataset_index];
	Random random = new Random();
	//    private HashMap<Integer, Double>[] kvs;
	HashMap<Integer, Double> kvs1;
	HashMap<Integer, Double>[] kvs;
	ArrayList<HashMap<ArrayKey, Integer>>[] frequencyMap = new ArrayList[2];
	ArrayList<HashMap<ArrayKey, Integer>>[] small_frequencyMap = new ArrayList[2];
	//	ArrayList<HashMap<ArrayKey, Integer>>[] frequencyMap_with_key = new ArrayList[2];
	ArrayList<HashMap<ArrayKey, Integer>>[] dist_1 = new ArrayList[2];
	//	ArrayList<HashMap<ArrayKey, Integer>>[] dist_before = new ArrayList[2];
	ArrayList<HashMap<ArrayKey, Integer>>[] frequency4mean = new ArrayList[2];
	ArrayList<HashMap<ArrayKey, Integer>>[] frequency4mean_short = new ArrayList[2];
	ArrayList<HashMap<ArrayKey, Integer>>[] frequency_2bit = new ArrayList[2];

	double[] empirical_eps_between = new double[iter];
	double[] empirical_eps_mean = new double[iter];
	//	Processor processor;
	int test_statistic = 0;
	int[] binNum = {4};

	int[][] perturb_st_predict = new int[2][binNum[0]];
	int[][] perturb_st_predict_mean = new int[2][binNum[0]];
	int n;
	int d;
	double user_id;
//	int[] threshold = {21, 21, 30, 74, 468, 4592, 1732, 1237, 556, 213, 2};
//	int[] threshold = {2187, 2355, 2868, 3934, 5876, 4225, 3073, 830, 147, 18, 0};

	int thre = 0;
	int[] mean_index = new int[iter];
	int[][] boolof1 = new int[2][iter];
	int[] binarysum = new int[iter];
	double[] empirical_eps_with_key = new double[9];
	int g_of_total = 8;
	double[][] outputdistribution = new double[2][g_of_total];
	double[][] outputdistribution_with_key = new double[2][g_of_total + 1];
	double[] empirical_eps_uehis_with_key = new double[iter];
	double[] empirical_eps = new double[g_of_total];
	double[] empirical_eps_uehis = new double[iter];
	//    int[] decimalValue = n ew int[iter];
	StringBuilder binaryString = new StringBuilder();
	int decimalValue;
	//
	int[] keyfre = {0, 0};
	int[][] hier = new int[2][binNum[0]];
	int bound = 2 ;
	int[][] hier_fre = new int[2][binNum[0]];

	HashMap<Integer, Double> kvs_ori = new HashMap<>();

	int[] perturb_encoding_tmp = new int[binNum[0]];

	double[] boundary = new double[binNum[0]];
	double[][] assignedMean = new double[2][binNum[0]  / 2];
	//    double[] boundary = {-1, -0.5, -0.5, 0.0, 0.0, 0.5, 0.5, 1.0};
	int[][] ptb = new int[2][iter];
	List<Map.Entry<Integer, Double>> kvlist = new ArrayList<>(kvs_ori.entrySet());
	double[] mean_whole = new double[]{0.0, 0.0, 0.0, 0.0};
	double[] count_whole = new double[binNum[0] / 2]; //{0.0, 0.0, 0.0, 0.0};
	public Histogram_KV_audit_compare() throws Exception{
		histogramTest();
//		}
	}


	public void histogramTest(){
//		int[] binNum = {4, 8, 12, 16};

//		double[] eps = {0.5, 1, 2, 4, 8, 16};
		double[] eps = {0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 4, 6};
		DecimalFormat df = new DecimalFormat("0.############");

		HashMap<Integer, Double>[] kvs = new HashMap[2];
		int method = 0;
		if(method==0)
			System.out.println("\n\t-------------  Histogram Test: PrivKVM*-GVPP: -------------");
		else if(method==1)
			System.out.println("\n\t-------------  Histogram Test: PrivKVM*-SW: -------------");
		else if(method==2)
			System.out.println("\n\t-------------  Histogram Test: PrivKVM*-PM: -------------");
		boundary[0] = -1;
		boundary[binNum[0]-1] = 1;
		double step = (double)2 / (binNum[0]/2);
		for(int i = 1; i < binNum[0] - 1; i = i + 2){
			boundary[i] = boundary[i - 1] + step;
			boundary[i + 1] = boundary[i];
		}
		int meanindex = 0;
		for(int i = 0; i < binNum[0] - 1; i = i + 2){
			assignedMean[0][meanindex] = (boundary[i] + boundary[i + 1])/2;
			assignedMean[1][meanindex] = (boundary[i] + boundary[i + 1])/2;
			meanindex++;
		}
//		assignedMean[0] = new double[]{-1.0, 0, 0,1.0};
//		assignedMean[1] = new double[]{-1.0, 0, 0,1.0};
		for(int i=0; i<binNum.length; i++){
			System.out.println("\n******************** bin_num = " + binNum[i] + " ********************");
			for(int j=0; j<eps.length; j++){
//				thre = threshold[j];
				double average_mse = 0;
				double[] histogram_per = null;
				int index = 0;
				kvs_ori.put(2, -1.0);

				kvs_ori.put(4, 1.0);
				meanindex = 0;
				for(int imean = 0; imean < binNum[0] - 1; imean = imean + 2){
					assignedMean[0][meanindex] = (boundary[imean] + boundary[imean + 1])/2;
					assignedMean[1][meanindex] = (boundary[imean] + boundary[imean + 1])/2;
					meanindex++;
				}
				kvlist.addAll(new HashMap<>(kvs_ori).entrySet());
				double[] histogram_real = new double[binNum[0] / 2];
//				{0.0, 0.0, 0.0, iter/2};
				histogram_real[binNum[0]/2-1] = iter/2;
				for(int ih = 0; ih < binNum[0] / 2 - 1; ih++){
					histogram_real[ih] = 0.0;
				}
				Map.Entry<Integer, Double> kv = kvlist.get(0);
				for(int iii = 0; iii < 2; iii++){
					frequencyMap[iii] = new ArrayList<>();
					dist_1[iii] = new ArrayList<>();
					frequency4mean[iii] = new ArrayList<>();
					frequencyMap[iii].add(new HashMap<>());
					dist_1[iii].add(new HashMap<>());
					frequency4mean[iii].add(new HashMap<>());
					frequency_2bit[iii] = new ArrayList<>();
					frequency_2bit[iii].add(new HashMap<>());
				}
				for(int tmp = 0; tmp < 5; tmp++)
					histogram_UE(0, binNum[0], eps[j], kv,-1);
//
				for(ArrayKey key : frequencyMap[0].get(0).keySet()) {
					if(!dist_1[0].get(0).containsKey(key))
						dist_1[0].get(0).put(key, frequencyMap[0].get(0).get(key));
					else{
						int value = dist_1[0].get(0).get(key);
						dist_1[0].get(0).put(key, frequencyMap[0].get(0).get(key) + value);
					}
				}
				kv = kvlist.get(1);
				for(int tmp = 0; tmp < 5; tmp++)
					histogram_UE(1, binNum[0], eps[j], kv,-1);
//
				for(ArrayKey key : frequencyMap[1].get(0).keySet()) {
					if(!dist_1[1].get(0).containsKey(key))
						dist_1[1].get(0).put(key, frequencyMap[1].get(0).get(key));
					else{
						int value = dist_1[1].get(0).get(key);
						dist_1[1].get(0).put(key, frequencyMap[1].get(0).get(key) + value);
					}
				}   //只是为了计算出所有数据都参与训练的时候数据的分布

				for(ArrayKey key : dist_1[1].get(0).keySet()) {
					dist_1[1].get(0).put(key, dist_1[1].get(0).get(key)/5);
				}
				for(ArrayKey key : dist_1[0].get(0).keySet()) {
					dist_1[0].get(0).put(key, dist_1[0].get(0).get(key)/5);
				}

				meanindex = 0;
				assignedMean = new double[2][binNum[0] / 2];
				for(int imean = 0; imean < binNum[0] - 1; imean = imean + 2) {
					assignedMean[0][meanindex] = (boundary[imean] + boundary[imean + 1]) / 2;
					assignedMean[1][meanindex] = (boundary[imean] + boundary[imean + 1]) / 2;
					meanindex++;
				}


				for(int user_id = 0; user_id < 5; user_id++){
					if(method==0) {
						histogram_per = histogram_PrivKVM(user_id, kvlist, eps[j], binNum[i], user_id);

					}
					double mse = 0;
					for(int l=0; l<binNum[i]/2; l++)
						mse += Math.pow(histogram_real[l]-histogram_per[l], 2)/(binNum[i]/2);

					average_mse += mse/d;

				}


			}

		}


	}



	public double[] get_histogram(int binNum, double[] value){
		double[] count = new double[binNum];
		double[] boundary = getBoundary(binNum);
		for(int i=0; i<value.length; i++){
			if(value[i] == Double.MAX_VALUE)
				continue;
			if(value[i]>=boundary[0] && value[i]<=boundary[1])  // in the first interval, the left boundary should be included
				count[0]++;
			for(int j=1; j<boundary.length-1; j++){
				if(value[i]>boundary[j] && value[i]<=boundary[j+1])
					count[j]++;
			}
		}
		return count;
	}

	// set the boundary with binNum bins, in ascending order
	public double[] getBoundary(int binNum){
		double[] boundary = new double[binNum+1];
		boundary[0] = -1;
		boundary[binNum] = 1;
		for(int i=1; i<binNum; i++)
			boundary[i] = boundary[0] + 2.0/binNum*i;
		return boundary;
	}




	public double PiecewiseMechanismPerturbation(double value, double epsilon){
		double res = 0;
		double C = (Math.exp(epsilon/2)+1)/(Math.exp(epsilon/2)-1);
		double p = (Math.exp(epsilon) - Math.exp(epsilon/2)) / (2*Math.exp(epsilon/2) + 2);
		double l = (C+1)/2*value - (C-1)/2;
		double r = l + C - 1;

		Random random = new Random();
		if(random.nextDouble()<p){   // sample a value from [l,r]
			res = getRandomNum(l,r);
		}else{
			res = getRandomNum(-1*C, C);
			while(res<=r && res>=l)
				res = getRandomNum(-1*C, C);
		}
		res /= C;
		return res;
	}

	public static double getRandomNum(double l, double r){
		Random random = new Random();
		double rnd = l + (r-l)*random.nextDouble();
		return rnd;
	}


	public double SquarWaveMechanismPerturbation(double value, double epsilon){
		double res = 0;

		// normalize a value in [-1, 1] to the range of [0, 1]
		value = (value+1)/2.0;
		double temp = Math.exp(epsilon);
		double b = (epsilon*temp-temp+1)/(2*temp*(temp-1-epsilon));
		double p = temp/(2*b*temp+1);
		double q = 1.0/(2*b*temp+1);

		Random random = new Random();
		if(random.nextDouble()<p){
			res = value - b + (2*b) * random.nextDouble();
		}else{
			res = -b + (2*b+1) * random.nextDouble();
			while(res>value-b && res<value+b)
				res = -b + (2*b+1) * random.nextDouble();
		}

		res = 2*(res+b)/(2*b+1)-1;

		return res;
	}


	public double[] getAllValuesAfterKeyPerturbation(int index, Map.Entry<Integer, Double> kvs, double epsilon, int binNum, int times){
		double p = Math.exp(epsilon)/(1.0+Math.exp(epsilon));
		if(times > 0)
			p=0.5;
		double[] cnt_v = new double[binNum];
		double v = kvs.getValue();


		double ratio_p = 0.5;

//		double[] value = new double[n];
		double[] value = new double[iter];
		int cnt = 0;
		for(int i = 0; i < iter; i++){
			int forkey =  0;
			if(times == 0) {
				if (index == 0) {
//					if (random.nextDouble() < ratio_p) {
					forkey = 1;
//					}
					if (forkey == 1) {
						if (random.nextDouble() < p) {
							keyfre[index] = keyfre[index] + 1;
						}
					}
				}
			}

			if(times == 0) {
				if (index == 1) {
					if(random.nextDouble() < ratio_p){
						forkey = 0;
					}
					if (forkey == 0) {
						if (random.nextDouble() > p) {
							keyfre[index] = keyfre[index] + 1;
						}
					}
				}
			}

			if(times == -1) {
				boolof1[index][i] = 1;
			}else if(times == 0){
				int randomValue =  0 ;
				if(random.nextDouble() < ratio_p){
					randomValue = 1;
				}
				boolof1[index][i] = randomValue;


			}

			if((boolof1[index][i]== 1)) {
				if (random.nextDouble() < p) {
					value[cnt++] = v;  // 1->1
					ptb[index][i] = 1;
				}
				else {
					value[cnt++] = Double.MAX_VALUE;  // 1->0
					//			}else{
					ptb[index][i] = 0;
				}
			}else{
				if (random.nextDouble() < p) {  // 0->0
					value[cnt++] = Double.MAX_VALUE;
					ptb[index][i] = 1;
				}
				else {
					int idx = random.nextInt(binNum/2);
					value[cnt++] = assignedMean[index][idx];
					ptb[index][i] = 1;
				}
			}
		}
		return value;

	}
	public double calculate_cv(ArrayList<HashMap<ArrayKey, Integer>> input_list, int threshold){
		int b_nun = 0;
		double sum_befor = 0.0;
		for (int num : input_list.get(0).values()) {
			if(num < threshold)
				continue;
			sum_befor += num;
			b_nun++;
		}
		sum_befor = sum_befor / b_nun;
		// Step 2: 计算方差
		double variance_sort = 0.0;
		for (int num : input_list.get(0).values()) {
			if(num < threshold)
				continue;
			variance_sort += Math.pow(num - sum_befor, 2);
		}
		variance_sort /= b_nun;
		variance_sort = Math.sqrt(variance_sort);
		double cov_sort = variance_sort / sum_befor;
		return cov_sort;
	}



	public double calculate_cv_ori(ArrayList<HashMap<ArrayKey, Integer>> input_list, int threshold){
		int b_nun = 0;
		double sum_befor = 0.0;
		for (int num : input_list.get(0).values()) {
			if(num < threshold)
				continue;
			sum_befor += num;
			b_nun++;
		}
		sum_befor = sum_befor / b_nun;

		double variance_sort = 0.0;
		for (int num : input_list.get(0).values()) {
			if (num < 0)
				continue;
			variance_sort += Math.pow(num - sum_befor, 2);
		}
		variance_sort /= b_nun;
//		variance_sort /= iter;
		variance_sort = Math.sqrt(variance_sort);
		double cov_sort = variance_sort / sum_befor;
		return cov_sort;
	}

	public double[] calculate_threshold(ArrayList<HashMap<ArrayKey, Integer>> change, double ori_cv, double eps, String str1){
		double cal_cv = 0.0;
		int num_thresh = 5;
		if(eps <= 0.2)
			cal_cv = 6;

		else if(eps <= 0.4)
//			cal_cv = 2.5;
			cal_cv = 6;
		else
			cal_cv = 6;

		double sum2 = 0.0;
		double variance2 = 0.0;
		int threshold = 0;
		TreeMultimap<Integer, String> s1 = TreeMultimap.create(
				Ordering.natural().reverse(), // 降序排列键
				Ordering.natural()
		);

		for (HashMap<ArrayKey, Integer> entry : change) {
			for(ArrayKey key : entry.keySet()) {
				s1.put(change.get(0).get(key), Arrays.toString(key.array));
			}
		}
		List<ArrayKey> intersection = new ArrayList<>();


		int num_tree = 0;
		double cov_tree1 = 0.0;
		int num_1 = 0;
		for(Map.Entry<Integer, String> entry : s1.entries()) {
			int num = entry.getKey();
			if(num == 1){
				num_1++;
			}
			if(num < num_thresh){
				break;
			}
			num_tree++;
			sum2 += num;
			sum2 /= num_tree;
			int t_tmp = 0;
			variance2 = 0.0;
			for (Map.Entry<Integer, String> entry_t : s1.entries()) {
				if (t_tmp < num_tree) {
					int num_t = entry_t.getKey();
					variance2 += Math.pow(num_t - sum2, 2);
					t_tmp++;
				} else {
					break;
				}
			}
			variance2 /= num_tree;
			variance2 = Math.sqrt(variance2);
			cov_tree1 = variance2 / sum2;
			if (cov_tree1 >= cal_cv) {
				num_tree--;
				break;
			} else {
				sum2 *= num_tree;
				String str = new String(entry.getValue());
				str = str.substring(1, str.length() - 1);
				String[] stringArray = str.split(",\\s*");
				int[] intArray = new int[stringArray.length];
				for (int i = 0; i < stringArray.length; i++) {
					intArray[i] = Integer.parseInt(stringArray[i]);
				}
				ArrayKey tmp_a = new ArrayKey(intArray);
				if(!intersection.contains(tmp_a)) {
					intersection.add(tmp_a);
				}
			}
		}
//		}

		int numtmp = 1;
		int flag = 0;
		for(Map.Entry<Integer, String> entry : s1.entries()) {
			threshold = entry.getKey();
			if(numtmp == num_tree) {
				break;
			}
			numtmp++;
			flag = 1;
		}

		numtmp = 0;
		double sum = 0.0;
		for (Map.Entry<Integer, String> entry : s1.entries()) {
			if(numtmp == num_tree)
				break;
			double value = entry.getKey();
			sum += value;
		}
		double mean = sum / num_tree;
		numtmp = 0;
		for (Map.Entry<Integer, String> entry : s1.entries()) {
			if(numtmp == num_tree)
				break;
			double value = entry.getKey();
			sum += Math.pow(value - mean, 2);
		}


		double[] result = new double[2];
		result[0] = threshold;
		result[1] = 0;

		return result;
	}
	public double[] histogram_PrivKVM(int index,  List<Map.Entry<Integer, Double>> kvs, double epsilon, int binNum, int times){
		double[] count_per = null;
		DecimalFormat df = new DecimalFormat("0.############");
		for(int i = 0; i < 2; i++){
			frequencyMap[i] = new ArrayList<>();
			small_frequencyMap[i] = new ArrayList<>();
			frequency4mean_short[i] = new ArrayList<>();
			frequencyMap[i].add(new HashMap<>());
			small_frequencyMap[i].add(new HashMap<>());
			frequency4mean_short[i].add(new HashMap<>());
			frequency_2bit[i] = new ArrayList<>();
			frequency_2bit[i].add(new HashMap<>());
		}
		keyfre = new int[2];
		hier_fre = new int[2][binNum];
		String uebeforformatted = new String();
		perturb_st_predict = new int[2][binNum];

		Map.Entry<Integer, Double> kv = kvs.get(0);

		count_per = histogram_UE(0, binNum, epsilon, kv, times);
		kv = kvs.get(1);
		count_per = histogram_UE(1, binNum, epsilon, kv, times);


		System.out.println("eps" + epsilon);
		empirical_eps_mean = new double[iter];
//		List<ArrayKey> intersection_small = new ArrayList<>();
		HashMap<ArrayKey, Integer> tmp = new HashMap<>();

		HashMap<ArrayKey, Integer> tmp2 = new HashMap<>();
		ArrayList<HashMap<ArrayKey, Integer>>[] inter_frequencyMap = new ArrayList[2];

		List<ArrayKey> calculate = new ArrayList<>();

		ArrayList<ArrayKey> intersection = new ArrayList<>();

		System.out.println(perturb_st_predict[0][test_statistic] );
		double[] bound1_pre = clopperPearson(perturb_st_predict[0][test_statistic] , iter, 0.05);
		double[] bound2_pre = clopperPearson(perturb_st_predict[1][test_statistic], iter, 0.05);
		double log1_pre = Math.abs(Math.log(bound1_pre[0] / bound2_pre[1]));
		double log2_pre = Math.abs(Math.log(bound1_pre[1] / bound2_pre[0]));
		System.out.println("ac的 ： " +  Math.max(log2_pre, log1_pre));
//		System.out.println("均值取值" +  assignedMean[0][0] + " " + assignedMean[0][binNum/2-1] + " " +  assignedMean[1][0] + " "+  assignedMean[1][binNum/2-1]);
		perturb_st_predict = new int[2][binNum];

		int id = 0;

		int min_add = 0;
		double etap = Math.exp(epsilon/2.0/5.0);
	double eta = 1.0 / (1 + etap);
		for (ArrayKey key : frequencyMap[0].get(0).keySet()) {
			int tmpp = 0;
			if (dist_1[0].get(0).get(key) == null) {
				tmpp = 0;
			} else {
				tmpp = dist_1[0].get(0).get(key);
			}
		frequency4mean[0].get(0).put(key, frequencyMap[0].get(0).get(key) - (int) Math.floor(eta * tmpp));
		}

		for (ArrayKey key : frequencyMap[1].get(0).keySet()) {
			int tmpp = 0;
			if (dist_1[1].get(0).get(key) == null) {
				tmpp = 0;
			} else {
				tmpp = dist_1[1].get(0).get(key);
			}
			frequency4mean[1].get(0).put(key, frequencyMap[1].get(0).get(key) - (int) Math.floor(eta * tmpp));
		}



		//todo 缩小的交集
		intersection = new ArrayList<>();
//		calculate = new ArrayList<>();
		inter_frequencyMap[0] = new ArrayList<>();
		inter_frequencyMap[1] = new ArrayList<>();
		inter_frequencyMap[0].add(new HashMap<>());
		inter_frequencyMap[1].add(new HashMap<>());

		tmp = small_frequencyMap[0].get(0);
		for (ArrayKey intmp : tmp.keySet()) {
			tmp2 = small_frequencyMap[1].get(0);
			//						for (int listi2 = 0; listi2 < frequencyMap[1].size(); listi2++) {
			//							Map<ArrayKey, Integer> tmp2 = frequencyMap[1].get(listi2);
			if (tmp2.containsKey(intmp)) {
				intersection.add(intmp);
				inter_frequencyMap[0].get(0).put(intmp, tmp.get(intmp));
				inter_frequencyMap[1].get(0).put(intmp, tmp2.get(intmp));
			}
		}

		double cv_cal1 = calculate_cv_ori(inter_frequencyMap[0], 0);
		double cv_cal2 = calculate_cv_ori(inter_frequencyMap[1], 0);

//		System.out.println("原来的cv： " + cv_cal1 + " " + cv_cal2);

		double[] ans1 = calculate_threshold(inter_frequencyMap[0], cv_cal1, epsilon , "normal");
		int thre_s_1 = (int)ans1[0];
		double[] ans2 = calculate_threshold(inter_frequencyMap[1], cv_cal2, epsilon, "normal");
		int thre_s_2 = (int)ans2[0];
//		System.out.println("z-score :"+ ans1[1] + " " + ans2[1]);
		double cv_s1 = calculate_cv(inter_frequencyMap[0], thre_s_1);
		double cv_s2 = calculate_cv(inter_frequencyMap[1], thre_s_2);

		empirical_eps_uehis = new double[iter];
		empirical_eps_between = new double[iter];



		double[] bound_before0 = new double[2];


		calculate = new ArrayList<>();

		intersection = new ArrayList<>();
//		calculate = new ArrayList<>();
		inter_frequencyMap[0] = new ArrayList<>();
		inter_frequencyMap[1] = new ArrayList<>();
		inter_frequencyMap[0].add(new HashMap<>());
		inter_frequencyMap[1].add(new HashMap<>());

		tmp = frequency4mean_short[0].get(0);
		for (ArrayKey intmp : tmp.keySet()) {
			tmp2 = frequency4mean_short[1].get(0);
			//						for (int listi2 = 0; listi2 < frequencyMap[1].size(); listi2++) {
			//							Map<ArrayKey, Integer> tmp2 = frequencyMap[1].get(listi2);
			if (tmp2.containsKey(intmp)) {
				intersection.add(intmp);
				inter_frequencyMap[0].get(0).put(intmp, tmp.get(intmp));
				inter_frequencyMap[1].get(0).put(intmp, tmp2.get(intmp));
			}
		}

		empirical_eps_uehis = new double[iter];
		empirical_eps_between = new double[iter];


		hier = new int[2][binNum];
		for(ArrayKey key : frequency4mean[0].get(0).keySet()){
			String arrayString = key.toString(); // 输出: [1, 2, 3, 4]
			String[] elements = arrayString.substring(1, arrayString.length() - 1).split(", ");

			// 转换为整数数组
			int[] array = Arrays.stream(elements).mapToInt(Integer::parseInt).toArray();


			for(int k = 0; k < array.length; k++) {
				hier[0][k] +=  (array[k] * frequency4mean[0].get(0).get(key));
			}
		}

		for(ArrayKey key : frequency4mean[1].get(0).keySet()){
			String arrayString = key.toString(); // 输出: [1, 2, 3, 4]
			String[] elements = arrayString.substring(1, arrayString.length() - 1).split(", ");

			// 转换为整数数组
			int[] array = Arrays.stream(elements).mapToInt(Integer::parseInt).toArray();


			for(int k = 0; k < array.length; k++) {
				hier[1][k] += (array[k] * frequency4mean[1].get(0).get(key));
			}
		}
		double[] eps_hier = new double[binNum];
//			for (int k = 0; k < binNum; k++) {
//				System.out.println(k + ": " + hier[0][k] + " " + hier[1][k]);
//			}
		for (int k = 0; k < binNum; k++) {
			if(hier[0][k] <= 0 || hier[1][k] <= 0)
				continue;
			double[] bound1 = clopperPearson(hier[0][k], iter, 0.05);
			double[] bound2 = clopperPearson(hier[1][k], iter, 0.05);
			double eps1 = Math.abs(Math.log(bound1[0] / bound2[1]));
			double eps2 = Math.abs(Math.log(bound1[1] / bound2[0]));
			double eps3 = Math.abs(Math.log((iter - bound1[1]) / (iter - bound2[0])));
			double eps4 = Math.abs(Math.log((iter - bound1[0]) / (iter - bound2[1])));
			eps_hier[k] = Math.max(eps2, eps1);
			eps_hier[k] = Math.max(eps_hier[k], eps3);
			eps_hier[k] = Math.max(eps_hier[k], eps4);
		}
		double max_eps_hier = Arrays.stream(eps_hier).max().getAsDouble();
		System.out.println("hier 4 mean: " + max_eps_hier);
		tmp = frequency4mean[0].get(0);
		inter_frequencyMap[0] = new ArrayList<>();
		inter_frequencyMap[1] = new ArrayList<>();
		inter_frequencyMap[0].add(new HashMap<>());
		inter_frequencyMap[1].add(new HashMap<>());
		intersection = new ArrayList<>();
		//todo 直接扰动完的inter
		for (ArrayKey intmp : tmp.keySet()) {
			tmp2 = frequency4mean[1].get(0);
			if (tmp2.containsKey(intmp)) {
				intersection.add(intmp);
				inter_frequencyMap[0].get(0).put(intmp, tmp.get(intmp));
				inter_frequencyMap[1].get(0).put(intmp, tmp2.get(intmp));
				//							}
			}
		}
		cv_cal1 = calculate_cv_ori(inter_frequencyMap[0], 0);
		cv_cal2 = calculate_cv_ori(inter_frequencyMap[1], 0);
//		System.out.println("原来的cv： " + cv_cal1 + " " + cv_cal2);

		ans1 = calculate_threshold(inter_frequencyMap[0], cv_cal1, epsilon, "mean");
		int thre_m1 = (int)ans1[0];
		ans2 = calculate_threshold(inter_frequencyMap[1], cv_cal2, epsilon, "mean");
		int thre_m2 = (int)ans2[0];

		cv_cal1 = calculate_cv(inter_frequencyMap[0], thre_m1);
		cv_cal2 = calculate_cv(inter_frequencyMap[1], thre_m2);


//		thre_m1 = 1;
//		thre_m2 = 1;
//		if(thre_m1 == 0 || thre_m2 == 0){
//			thre_m2 = 1;
//			thre_m1 = 1;
//		}
		double cv_m1 = calculate_cv(inter_frequencyMap[0], thre_m1);
		double cv_m2 = calculate_cv(inter_frequencyMap[1], thre_m2);
		double log1 = 0.0, log2 = 0.0;
//		int size1 = frequency4mean[0].get(0).size();
//		System.out.println(intersection.size());
		for (ArrayKey key : intersection) {

			double[] bound_before1 = new double[2];
			double[] mean1 = new double[2];
			double[] mean2 = new double[2];

			double log_m2 = 0.0, log_m1 = 0.0;
//			System.out.println(thre_m2 + " " + thre_m1);
			if (frequency4mean[1].get(0).get(key) <  thre_m2|| (frequency4mean[0].get(0).get(key) < thre_m1)) {
//			if (frequency4mean[1].get(0).get(key) <  5|| (frequency4mean[0].get(0).get(key) < 5)) {
//				empirical_eps_mean[id] = 0.0;
				continue;
			}else{
//				System.out.println(frequency4mean[1].get(0).get(key));
				mean1 = clopperPearson((int) frequency4mean[1].get(0).get(key), iter, 0.05);
//				System.out.println(frequency4mean[0].get(0).get(key));
				mean2 = clopperPearson((int) frequency4mean[0].get(0).get(key), iter, 0.05);
				log_m1 = Math.abs(Math.log(mean1[0] / mean2[1]));
				log_m2 = Math.abs(Math.log(mean1[1] / mean2[0]));
//				if(Double.isFinite(log_m1) || Double.isFinite(log_m2)){
//					continue;
//				}


			}
			empirical_eps_mean[id] = Math.abs(Math.max(log_m1, log_m2));
//			if(empirical_eps_mean[id] > epsilon){
////				id = id + 1;
//				continue;
//			}
			id = id + 1;
			calculate.add(key);
		}

		id = 0;


		uebeforformatted = df.format(Arrays.stream(empirical_eps_mean).max().getAsDouble());
		OptionalInt maxIndex = IntStream.range(0, empirical_eps_mean.length)
				.reduce((i, j) -> empirical_eps_mean[i] > empirical_eps_mean[j] ? i : j);

		int index_max = maxIndex.getAsInt();
		int tmp_index = 0;
		ArrayKey maxkey = new ArrayKey(perturb_encoding_tmp);
		for(ArrayKey key : calculate){
			maxkey = key;

			if(tmp_index== index_max){
				break;
			}
			tmp_index++;
		}
		System.out.println("mean " + uebeforformatted + " " +thre_m1 + " " + thre_m2+ " " + cv_cal1 + " " + cv_cal2 + " " + cv_m1 + " " + cv_m2 + " " +frequency4mean[1].get(0).get(maxkey) + " " + frequency4mean[0].get(0).get(maxkey));

		if(times == 0) {
			double[] bound_k1 = clopperPearson(keyfre[0], iter, 0.05);
			double[] bound_k2 = clopperPearson(keyfre[1], iter, 0.05);
			double epsk1 = Math.abs(Math.log(bound_k1[0] / bound_k2[1]));
			double epsk2 = Math.abs(Math.log(bound_k1[1] / bound_k2[0]));
			double epskey = Math.max(epsk1, epsk2);
//		epskey = Arrays.stream(epskey).max().getAsDouble();
			System.out.println("key : " + epskey);
		}
		eps_hier = new double[binNum];
//		for (int k = 0; k < binNum; k++) {
//			System.out.println((double)hier[0][k]/iter);
//		}
		for (int k = 0; k < binNum; k++) {
			double[] bound1 = clopperPearson(hier_fre[0][k], iter, 0.05);
			double[] bound2 = clopperPearson(hier_fre[1][k], iter, 0.05);
			double eps1 = Math.abs(Math.log(bound1[0] / bound2[1]));
			double eps2 = Math.abs(Math.log(bound1[1] / bound2[0]));
			double eps3 = Math.abs(Math.log((iter - bound1[1]) / (iter - bound2[0])));
			double eps4 = Math.abs(Math.log((iter - bound1[0]) / (iter - bound2[1])));
			eps_hier[k] = Math.max(eps2, eps1);
			eps_hier[k] = Math.max(eps_hier[k], eps3);
			eps_hier[k] = Math.max(eps_hier[k], eps4);
		}
		max_eps_hier = Arrays.stream(eps_hier).max().getAsDouble();
		System.out.println("hier : " + max_eps_hier);
		small_frequencyMap[0] = new ArrayList<>();
		small_frequencyMap[1] = new ArrayList<>();
		small_frequencyMap[0].add(new HashMap<>());
		small_frequencyMap[1].add(new HashMap<>());




		//正常的

		int thre1 = 0, thre2 = 0;
		tmp = frequencyMap[0].get(0);

		inter_frequencyMap[0] = new ArrayList<>();
		inter_frequencyMap[1] = new ArrayList<>();
		inter_frequencyMap[0].add(new HashMap<>());
		inter_frequencyMap[1].add(new HashMap<>());
		intersection = new ArrayList<>();
		//todo 直接扰动完的inter
		for (ArrayKey intmp : tmp.keySet()) {
			tmp2 = frequencyMap[1].get(0);
			if (tmp2.containsKey(intmp)) {
				intersection.add(intmp);
				inter_frequencyMap[0].get(0).put(intmp, tmp.get(intmp));
				inter_frequencyMap[1].get(0).put(intmp, tmp2.get(intmp));
				//							}
			}
		}
		cv_cal1 = calculate_cv_ori(inter_frequencyMap[0], 0);
		cv_cal2 = calculate_cv_ori(inter_frequencyMap[1], 0);
//		System.out.println("ue formatted 原来的cv： " + cv_cal1 + " " + cv_cal2 );

		int[] tmp_his = new int[iter];
		int ind = 0;
		for(ArrayKey thm: frequencyMap[0].get(0).keySet()){
			int inte = frequencyMap[0].get(0).get(thm);
			tmp_his[ind] = inte;
			ind++;
		}



		ans1 = calculate_threshold(inter_frequencyMap[0], cv_cal1, epsilon, "normal");
		thre1 = (int)ans1[0];
//		thre1 = 1;
		ans2 = calculate_threshold(inter_frequencyMap[1], cv_cal2, epsilon, "normal");
		thre2 = (int)ans2[0];
//		thre2 = 1;
		double cv_1 = calculate_cv(inter_frequencyMap[0], thre1);
		double cv_2 = calculate_cv(inter_frequencyMap[1], thre2);

//		ArrayList<ArrayKey>intersection1 = new ArrayList<>();
//		ArrayList<ArrayKey>intersection2 = new ArrayList<>(); #todo 其实不需要这个 我只要计算出他的那个阈值就可以
		empirical_eps_uehis = new double[iter];
		calculate = new ArrayList<>();

		for (ArrayKey key : intersection) {

			double[] bound_before1 = new double[2];
			double[] mean1 = new double[2];
			double[] mean2 = new double[2];

			double[] bound1 = new double[2];
			double[] bound2 = new double[2];
			if (frequencyMap[1].get(0).get(key)  < thre2|| frequencyMap[0].get(0).get(key) < thre1) {
				continue;

			} else {
				bound1 = clopperPearson((int) frequencyMap[1].get(0).get(key), iter, 0.05);
				bound2 = clopperPearson((int) frequencyMap[0].get(0).get(key), iter, 0.05);
			}
			log1 = Math.abs(Math.log(bound1[0] / bound2[1]));
			log2 = Math.abs(Math.log(bound1[1] / bound2[0]));
			empirical_eps_uehis[id] = Math.max(log1, log2);
			empirical_eps_uehis[id] = Math.abs(empirical_eps_uehis[id]);



			calculate.add(key);
			id = id + 1;

		}

		id = 0;

		uebeforformatted = df.format(Arrays.stream(empirical_eps_uehis).max().getAsDouble());
		maxIndex = IntStream.range(0, empirical_eps_uehis.length)
				.reduce((i, j) -> empirical_eps_uehis[i] > empirical_eps_uehis[j] ? i : j);

		index_max = maxIndex.getAsInt();
		tmp_index = 0;
		maxkey = new ArrayKey(perturb_encoding_tmp);
		for(ArrayKey key : calculate){
			maxkey = key;

			if(tmp_index== index_max){
				break;
			}
			tmp_index++;
		}

		System.out.println("ue formatted " + uebeforformatted + " " +cv_cal1 + " " + cv_cal2+ " " + cv_1 + " " + cv_2 + " " +frequencyMap[1].get(0).get(maxkey) + " " + frequencyMap[0].get(0).get(maxkey));

		id = 0;
		empirical_eps_uehis = new double[iter];
		for (ArrayKey key : frequency_2bit[0].get(0).keySet()) {

			double[] bound1 = new double[2];
			double[] bound2 = new double[2];

			bound1 = clopperPearson((int) frequency_2bit[1].get(0).get(key), iter, 0.05);
			bound2 = clopperPearson((int) frequency_2bit[0].get(0).get(key), iter, 0.05);

			log1 = Math.abs(Math.log(bound1[0] / bound2[1]));
			log2 = Math.abs(Math.log(bound1[1] / bound2[0]));
			empirical_eps_uehis[id] = Math.max(log1, log2);
			empirical_eps_uehis[id] = Math.abs(empirical_eps_uehis[id]);


			calculate.add(key);
			id = id + 1;

		}
		uebeforformatted = df.format(Arrays.stream(empirical_eps_uehis).max().getAsDouble());
		maxIndex = IntStream.range(0, empirical_eps_uehis.length)
				.reduce((i, j) -> empirical_eps_uehis[i] > empirical_eps_uehis[j] ? i : j);

		index_max = maxIndex.getAsInt();
		tmp_index = 0;
		maxkey = new ArrayKey(perturb_encoding_tmp);
		for(ArrayKey key : calculate){
			maxkey = key;

			if(tmp_index== index_max){
				break;
			}
			tmp_index++;
		}

		System.out.println("ue 2 bit " + uebeforformatted); // + " " +cv_cal1 + " " + cv_cal2+ " " + cv_1 + " " + cv_2 + " " +frequencyMap[1].get(0).get(maxkey) + " " + frequencyMap[0].get(0).get(maxkey));


		//todo 缩小的交集
		intersection = new ArrayList<>();


		return count_per;
	}


	public double[] histogram_UE(int index, int binNum, double epsilon,  Map.Entry<Integer, Double> kvs, int times){

		double[] values = getAllValuesAfterKeyPerturbation(index, kvs, epsilon/2.0, binNum, times); //todo 是不是先把数据给扰动了啊？


		double[] cnt_v = new double[binNum];
		double[] c_m = new double[binNum];

		//在这里就可以加入对values分布的统计了吧 不行 还得看那个扰动之后的unary的情况 但是key的话要考虑加入进来
		double[] noisyCount = new double[binNum];
//        double[] noisy_Count_fuben = new double[binNum];
		double[] noisyCount_with_key = new double[binNum + 1];
		int[] perturb_encoding_with_key = new int[binNum + 1];

		for(int i=0; i<values.length; i++) {
			int[] perturb_encoding = valuePert_UnaryEncoding_histogram_4mean(binNum, epsilon/2.0/5, values[i], i);
			for (int k = 0; k < noisyCount.length; k++)
				noisyCount[k] += perturb_encoding[k];  //这个就是没有扰校正的那个OUE的情况
			ArrayList<Integer> perturb_st_predict_idx = new ArrayList<>();
			for(int k = 0; k < noisyCount.length; k++){
				if(perturb_encoding[k] == 1)
					perturb_st_predict_idx.add(k);
			}

			int predict_idx = 0;
			if(perturb_st_predict_idx.size() == 0){
				predict_idx = random.nextInt(noisyCount.length);
			}else {
				predict_idx = perturb_st_predict_idx.get(random.nextInt(perturb_st_predict_idx.size()));
			}
			perturb_st_predict[index][predict_idx] += 1;
			// todo：从是1的里面随机抽取一个index进行保存

			ArrayKey key = new ArrayKey(perturb_encoding);
//			ArrayKey key = new ArrayKey(mix_perturb_encoding);
			if (!frequencyMap[index].get(0).containsKey(key)) {
				frequencyMap[index].get(0).put(key, 1);
			} else {
				int value = frequencyMap[index].get(0).get(key);
				frequencyMap[index].get(0).put(key, value + 1);
			}
			int[] sele = {perturb_encoding[0], perturb_encoding[binNum - 1]};
			ArrayKey select = new ArrayKey(sele);
			if (!frequency_2bit[index].get(0).containsKey(select)) {
				frequency_2bit[index].get(0).put(select, 1);
			} else {
				int value = frequency_2bit[index].get(0).get(select);
				frequency_2bit[index].get(0).put(select, value + 1);
			}

			for(int k = 0; k < hier_fre[index].length; k++){
				hier_fre[index][k] += perturb_encoding[k];
			}

		}


		int numberOfValues = 0;
		for(double v : values)
			if(v<Double.MAX_VALUE)
				numberOfValues++; //统计每一个value它的个数
		calibrate(noisyCount, 2, numberOfValues, epsilon/2.0/5);  //这个是单独校正count todo：和calibrate_oue相同的  // epsilon for Unary Encoding with sensitivity of 2

		int j = 0;
		for(int i = 0; i < noisyCount.length ; i = i+ 2){
			assignedMean[index][j] = (boundary[i+1] * noisyCount[i+1] + boundary[i]*noisyCount[i]) / (noisyCount[i] + noisyCount[i + 1]) ;
			j = j + 1;
		}
		//todo：确实是聚合部分的校正 是每个bucket的校正
		int counti = 0;
		for(int i=0; i<noisyCount.length; i = i + 2) { //第二次校正 应该是聚合部分的校正吧
//			|
			count_whole[counti] = rangeCountCalibration(noisyCount[i] + noisyCount[i + 1], 0, epsilon/2.0/5, binNum / 2 + 1, numberOfValues);
			counti = counti + 1;
		}
		return cnt_v; //count_whole;  //这个是用公式17进行校正的
	}



	//todo 但是你这个计算你也只是一个数据去计算啊 你直方图去计算的话只有一个左边界的就算是 去进行扰动的
	// due to key sampling
	public double rangeCountCalibration(double v, int key, double epsilon, int g, int total){
//		double frequency = this.processor.frequency[key];
		double p = Math.exp(epsilon)/(1+Math.exp(epsilon));
		double frequency = (double)total / iter;
		frequency = (p - 1 + frequency) / (2 * p - 1);
		double res = iter*v/p/total - iter*(1-frequency)*(1-p)/(g-1)/p;
		if(res<0)
			res = 0;
		if(res>total)
			res = total;
		return res;
	}


	public int valuePert_GRR_histogram(int binNum, double ep, double value){
		if(value == Double.MAX_VALUE)
			return Integer.MAX_VALUE;
		int idx_bin = rangeDiscretization_histogram(value, binNum);  //看这个值是落到了哪个bucket中的
		return GRR(idx_bin, binNum, ep);  //把这个idx给它进行扰动的
	}


	public static int GRR(int index, int num, double epsilon){
		double exp = Math.exp(epsilon);
		double p = exp/(exp+num-1);
		Random random = new Random();
		if(random.nextDouble() < p)
			return index;
		int idx = random.nextInt(num);
		while(idx==index)
			idx = random.nextInt(num); //哦还可以这样找不同的啊？
		return idx;
	}


	public int rangeDiscretization_histogram(double value, int binNum){
		double[] boundary = getBoundary(binNum);
		if(value>=boundary[0] && value<=boundary[1])
			return 0;
		int index = 1;
		for(int i=1; i<binNum; i++){
			if(value>boundary[i] && value<=boundary[i+1]){
				index = i;
				break;
			}
		}
		return index;
	}

	public int rangeDiscretization_histogram_4mean(double value, int binNum, int in){
		double p = 0.0;
		Random random = new Random();
		if(value>=boundary[0] && value<=boundary[1]){
			p = (boundary[1] - value) / (boundary[1] - boundary[0]);
			if(random.nextDouble() < p) {
				mean_index[in] = 0;
				return 0;

			}else{
				mean_index[in] = 1;
				return 1;
			}
		}
		int index = 1;
		for(int i=2; i<binNum; i = i + 2){
			if(value>boundary[i] && value<=boundary[i+1]){
				p = (boundary[i + 1] - value) / (boundary[i+1] - boundary[i]);
				if(random.nextDouble() < p){
//                    index = 2 * i;
					index = i;
				}else{
//                    index = 2 * i + 1;
					index = i+1;
				}
				break;
			}
		}
		mean_index[in] = index;
		return index;
	}


	public static double[] valueCountCalibration(int[] perturb_index, int num, double ep){
		double[] count = new double[num];
		int cnt = 0;
		for(int i=0; i<perturb_index.length; i++){
			if(perturb_index[i]<Integer.MAX_VALUE){
				count[perturb_index[i]]++; //每个perturb的index的个数
				cnt++;  //总的样本个数
			}
		}
		calibrate(count, num, cnt, ep);
		return count;
	}

	public static void calibrate(double[] count, int k, int N, double ep){
		double exp = Math.exp(ep);
		double p = exp / (k-1+exp);
		for(int i=0; i<count.length; i++){
			count[i] = (int)Math.round( N*(p-1)/(k*p-1.0) + count[i]*(k-1.0)/(k*p-1.0) );  //todo：这个公式真的不知道是哪里来的啊
//			double tmp = (2 *(exp + 1)* count[i] - 2 * N ) / (exp - 1);
			if(count[i]>N)
				count[i] = N;
			if(count[i]<0)
				count[i] = 0;
		}
	}


	public int[] valuePert_UnaryEncoding_histogram_4mean(int binNum, double ep, double value, int in){
		if(value == Double.MAX_VALUE)
			return new int[binNum];
		int idx = rangeDiscretization_histogram_4mean(value, binNum, in); // 1->left, 2->right, 0->others  //让他转换成直方图的形式
		return UnaryEncoding_4mean(idx, binNum , ep);  //返回的是编码
	}

	public int[] valuePert_UnaryEncoding_histogram(int binNum, double ep, double value){
		if(value == Double.MAX_VALUE)
			return new int[binNum];
		int idx = rangeDiscretization_histogram(value, binNum); // 1->left, 2->right, 0->others  //让他转换成直方图的形式
		return UnaryEncoding(idx, binNum, ep);  //返回的是编码
	}
	public static double[] clopperPearson(int k, int n, double alpha) {
		double lowerBound, upperBound;


		BetaDistribution betaLower = new BetaDistribution(k, n - k + 1);
		lowerBound = betaLower.inverseCumulativeProbability(alpha / 2);


		BetaDistribution betaUpper = new BetaDistribution(k + 1, n - k);
		upperBound = betaUpper.inverseCumulativeProbability(1 - alpha / 2);


		return new double[]{lowerBound, upperBound};
	}

	public static int[] UnaryEncoding_4mean(int index, int num, double ep){
		double exp = Math.exp(ep);
		double p = exp/(exp+1);
		int[] encoding = new int[num];
		Random random = new Random();
//        int index_4mean = index / 2;
		for(int i=0; i<num; i++){
			if(i==index){
				if(random.nextDouble()<=p) //生成一个随机值 (0,1]范围之内的
					encoding[i] = 1;
			}else{
				if(random.nextDouble()>p)
					encoding[i] = 1;
			}
		}
		return encoding;
	}

	public static int[] UnaryEncoding(int index, int num, double ep){
		double exp = Math.exp(ep);
		double p = exp/(exp+1);
		int[] encoding = new int[num];
		Random random = new Random();
		for(int i=0; i<num; i++){
			if(i==index){
				if(random.nextDouble()<=p) //生成一个随机值 (0,1]范围之内的
					encoding[i] = 1;
			}else{
				if(random.nextDouble()>p)
					encoding[i] = 1;
			}
		}
		return encoding;
	}

	public static void main(String[] args) throws Exception{
		new Histogram_KV_audit_compare();
	}
}

class ArrayKey {
	int[] array;

	public ArrayKey(int[] array) {
//        if (array.length != 3) {
//            throw new IllegalArgumentException("Array must be of length 4");
//        }
		this.array = array.clone(); // 避免外部修改数组
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj) return true;
		if (obj == null || getClass() != obj.getClass()) return false;
		ArrayKey other = (ArrayKey) obj;
		return Arrays.equals(array, other.array);
	}

	@Override
	public int hashCode() {
		return Arrays.hashCode(array);
	}

	@Override
	public String toString() {
		return Arrays.toString(array);
	}
}

