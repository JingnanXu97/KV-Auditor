package th_w;

//import kv_whole.ArrayKey;
import org.apache.commons.math3.distribution.BetaDistribution;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Random;
import java.util.List;
import java.util.*;
import java.util.stream.Collectors;

public class oue {


    int iter = 100000000;
    int v1 = 1;
    int v2 = 2;
    int binNum = 8;
    ArrayList<HashMap<List<Double>, Integer>>[] frequencyMap = new ArrayList[2];
    ArrayList<HashMap<List<Double>, Integer>>[] frequency_2bit = new ArrayList[2];

    double[] epsilon = {0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 4, 6};

    public oue() throws Exception{
        audit();
    }


    public int getRandomChoice(double[] ue_val, int k) {
        // Sum the elements of ue_val
        int sum = 0;
        for (double val : ue_val) {
            sum += val;
        }

        Random random = new Random();

        // If the sum is 0, return a random integer between 0 and k-1
        if (sum == 0) {
            return random.nextInt(k);
        } else {
            // Create a list to store indices where ue_val[i] == 1
            ArrayList<Integer> indices = new ArrayList<>();
            for (int i = 0; i < ue_val.length; i++) {
                if (ue_val[i] == 1) {
                    indices.add(i);
                }
            }

            // Choose a random index from the collected indices
            return indices.get(random.nextInt(indices.size()));
        }
    }



    public void audit(){
        for(int eps = 0; eps < epsilon.length; eps++) {
            frequencyMap[0] = new ArrayList<>();
            frequencyMap[0].add(new HashMap<>());
            frequencyMap[1] = new ArrayList<>();
            frequencyMap[1].add(new HashMap<>());
            frequency_2bit[0] = new ArrayList<>();
            frequency_2bit[0].add(new HashMap<>());
            frequency_2bit[1] = new ArrayList<>();
            frequency_2bit[1].add(new HashMap<>());
            double[] pert1 = new double[binNum];
            double[] pert2 = new double[binNum];
            double[] return1 = new double[binNum];
            double[] return2 = new double[binNum];
            double[] emp_eps = new double[binNum];
            double[] emp_int = new double[iter];
            int count1 = 0;
            int count2 = 0;
            int test_statistic = 1;
            for (int i = 0; i < iter; i++) {
                return1 = OUE(v1,  epsilon[eps]);
                List<Double> listKey = Arrays.stream(return1)
                        .boxed() 
                        .collect(Collectors.toList());

//                double[] tmp = {listKey.get(0), listKey.get(listKey.size()-1)};
                double[] tmp = {listKey.get(1), listKey.get(2)};
                List<Double> tmplk = Arrays.stream(tmp)
                        .boxed() 
                        .collect(Collectors.toList());
                if (!frequency_2bit[0].get(0).containsKey(tmplk)) {
                    frequency_2bit[0].get(0).put(tmplk, 1);
                }else{
                    int value = frequency_2bit[0].get(0).get(tmplk);
                    frequency_2bit[0].get(0).put(tmplk, value + 1);
                }

                 if (!frequencyMap[0].get(0).containsKey(listKey)) {
                    frequencyMap[0].get(0).put(listKey, 1);
                }else{
                    int value = frequencyMap[0].get(0).get(listKey);
                    frequencyMap[0].get(0).put(listKey, value + 1);
                }
              for(int j = 0; j < binNum; j++) {
                  pert1[j] += return1[j];
              }
                return2 = OUE(v2,  epsilon[eps]);
                listKey = Arrays.stream(return2)
                        .boxed() 
                        .collect(Collectors.toList());
              if (!frequencyMap[1].get(0).containsKey(listKey)) {
                    frequencyMap[1].get(0).put(listKey, 1);
                }else{
                    int value = frequencyMap[1].get(0).get(listKey);
                    frequencyMap[1].get(0).put(listKey, value + 1);
                }
             double[] ttmp ={listKey.get(1), listKey.get(2)};
               tmplk = Arrays.stream(ttmp)
                        .boxed()
                        .collect(Collectors.toList());
                if (!frequency_2bit[1].get(0).containsKey(tmplk)) {
                    frequency_2bit[1].get(0).put(tmplk, 1);
                }else{
                    int value = frequency_2bit[1].get(0).get(tmplk);
                    frequency_2bit[1].get(0).put(tmplk, value + 1);
                }
                for(int j = 0; j < binNum; j++) {
                    pert2[j] += return2[j];
                }
            }
            int id = 0;
            Map<List<Double>, Integer> intersection = new HashMap<>();
            HashMap<List<Double>, Integer> map1 = frequencyMap[0].get(0);
            HashMap<List<Double>, Integer> map2 = frequencyMap[1].get(0);

            for (Map.Entry<List<Double>, Integer> entry : map1.entrySet()) {
                List<Double> key = entry.getKey();
                if (map2.containsKey(key)) {
                    intersection.put(key, Math.min(entry.getValue(), map2.get(key)));
                }
            }




            for (Map.Entry<List<Double>, Integer> entry : intersection.entrySet()) {

                List<Double> listKey = entry.getKey();
                if(frequencyMap[0].get(0).get(listKey) < 1 || frequencyMap[1].get(0).get(listKey) < 1) {
                    continue;
                }
                double[] cp1 = clopperPearson((int)frequencyMap[0].get(0).get(listKey), iter, 0.05);
                double[] cp2 = clopperPearson((int)frequencyMap[1].get(0).get(listKey), iter, 0.05);
                double log_m1 = Math.abs(Math.log(cp1[0] / cp2[1]));
                double log_m2 = Math.abs(Math.log(cp1[1] / cp2[0]));
                emp_int[id] = Math.abs(Math.max(log_m1, log_m2));
                id = id +1;
            }
            System.out.println("hori: " + Arrays.stream(emp_int).max().getAsDouble());
            emp_int = new double[iter];

            intersection = new HashMap<>();
            map1 = frequency_2bit[0].get(0);
            map2 = frequency_2bit[1].get(0);
            for (Map.Entry<List<Double>, Integer> entry : map1.entrySet()) {
                List<Double> key = entry.getKey();
                if (map2.containsKey(key)) {
                    intersection.put(key, Math.min(entry.getValue(), map2.get(key)));
                }
            }

            id = 0;
            for (Map.Entry<List<Double>, Integer> entry : map1.entrySet()) {

                List<Double> listKey = entry.getKey();
                if(frequency_2bit[0].get(0).get(listKey) < 1 || frequency_2bit[1].get(0).get(listKey) < 1) {
                    continue;
                }
                double[] cp1 = clopperPearson((int)frequency_2bit[0].get(0).get(listKey), iter, 0.05);
                double[] cp2 = clopperPearson((int)frequency_2bit[1].get(0).get(listKey), iter, 0.05);
                double log_m1 = Math.abs(Math.log(cp1[0] / cp2[1]));
                double log_m2 = Math.abs(Math.log(cp1[1] / cp2[0]));
                emp_int[id] = Math.abs(Math.max(log_m1, log_m2));
                id = id +1;
            }
            System.out.println("2 key: " + Arrays.stream(emp_int).max().getAsDouble());
        }
    }


    public double[] OUE(int value, double epsilon) {
        double[] ue = new double[binNum];
        for(int i = 0; i < binNum; i++){
            if(value == i)
                ue[i] = 1;
        }
        double[] noiseue = new double[binNum];
        double p = 1.0/2;
        double q = 1 / (Math.exp(epsilon) + 1);
        for(int i = 0; i < binNum; i++){
            Random random = new Random();
            if(ue[i] != 1){
                double tmp = random.nextDouble();
                if(tmp <= p)
                    noiseue[i] = 1;
            }else{
                double tmp = random.nextDouble();
                if(tmp <= q)
                    noiseue[i] = 1;
            }
        }

        return noiseue;
    }

    public static double[] clopperPearson(int k, int n, double alpha) {
        double lowerBound, upperBound;


        BetaDistribution betaLower = new BetaDistribution(k, n - k + 1);
        lowerBound = betaLower.inverseCumulativeProbability(alpha / 2);


        BetaDistribution betaUpper = new BetaDistribution(k + 1, n - k);
        upperBound = betaUpper.inverseCumulativeProbability(1 - alpha / 2);
        return new double[]{lowerBound, upperBound};
    }
    public static void main(String[] args)  throws Exception {
        new oue();

    }

}

