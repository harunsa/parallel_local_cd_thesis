#include <gtest/gtest.h>
#include "adj_list.h"


TEST(FileReaderTest, smallFile){
    //test size
    AdjList graph;
    graph.addFromFile("../../google_test/UnitTestInputs/input1");
    ASSERT_EQ(graph.getSize(), 7);
}

TEST(EdgeInsertion, smallFile){
    AdjList graph;
    graph.addFromFile("../../google_test/UnitTestInputs/input1");
    // insert 1 edge, test if its inserted
    EXPECT_EQ(graph.getInnerTblCount(6), 2);
    // insert 1 edge, test if its undirected insertion
    EXPECT_EQ(graph.getEdgeCount(6), 2);
    // insert 2 identical Edges
    EXPECT_EQ(graph.getInnerTblCount(1), 2);
    EXPECT_EQ(graph.getEdgeCount(1), 2);
    // 1 different Edge at same timestamp
    EXPECT_EQ(graph.getInnerTblCount(2),4);
    EXPECT_EQ(graph.getEdgeCount(2),4);
    //insert 1 different Edge at different Timestamp
    EXPECT_EQ(graph.getInnerTblCount(4),2);
    EXPECT_EQ(graph.getEdgeCount(4),2);
    EXPECT_EQ(graph.getInnerTblCount(5),2);
    EXPECT_EQ(graph.getEdgeCount(5),2);
    // insert 2 different Edges with same source(1), at same timestamp
    EXPECT_EQ(graph.getInnerTblCount(7),3);
    EXPECT_EQ(graph.getDestSize(7,1),2);
    EXPECT_EQ(graph.getEdgeCount(7), 4);
    // insert 2 different Edges with same destination(2), at same timestamp
    EXPECT_EQ(graph.getInnerTblCount(8), 3);
    EXPECT_EQ(graph.getDestSize(8,2), 2);
    EXPECT_EQ(graph.getEdgeCount(8), 4);

}

TEST(EdgeDeletion, smallFile){
    AdjList graph;
    // delete Edge, same timestamp, inOrder
    graph.addFromFile("../../google_test/UnitTestInputs/input2");
    EXPECT_EQ(graph.getEdgeCount(1), 0);
    //delete Edge, same timestamp, outOfOrder
    EXPECT_EQ(graph.getEdgeCount(2), 0);
    // Edge to be deleted doesnt exist
    EXPECT_EQ(graph.getEdgeCount(3),2);
    // Edge to be deleted doesnt exist, duplicate deletes
    EXPECT_EQ(graph.getEdgeCount(4),0);
    // Edge to be deleted does exist, duplicate deletes
    EXPECT_EQ(graph.getEdgeCount(5),0);
    //delete Edge at source
    EXPECT_EQ(graph.getEdgeCount(6),2);
    EXPECT_EQ(graph.getDestSize(6,1),1);
    EXPECT_EQ(graph.getInnerTblCount(6),2);
    //delete edge at destination
    EXPECT_EQ(graph.getEdgeCount(7),2);
    EXPECT_EQ(graph.getDestSize(7,1),1);
    EXPECT_EQ(graph.getInnerTblCount(7),2);
}

//test findEdge
TEST(EdgeFind, smallFile){
    AdjList graph;
    graph.addFromFile("../../google_test/UnitTestInputs/input1");
    //find edge, it exists
    EXPECT_TRUE(graph.findEdge(1,2,1));
    //find edge, search for a non-existend Edge
    EXPECT_FALSE(graph.findEdge(1,3,3));
    //edge exists but in diff timestamp, out of range
    EXPECT_FALSE(graph.findEdge(3,2,5,7));
    //in range but edge does not exist
    EXPECT_FALSE(graph.findEdge(4,6,1,4));
    //in range and exists
    EXPECT_TRUE(graph.findEdge(1,2,0,10));
}


TEST(TimeConsumption,smallFile){
    //larger file
    AdjList graph;
    auto t1 = std::chrono::high_resolution_clock::now();
    graph.addFromFile("../../google_test/UnitTestInputs/input1");
    auto t2 = std::chrono::high_resolution_clock::now();
    auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    EXPECT_LT(ms_int.count(),500);
}

TEST(RangeQuery, smallFile){
    AdjList graph;
    graph.addFromFile("../../google_test/UnitTestInputs/input1");

    //out of range
    graph.rangeQuery(9,15, [](uint64_t timestamp,uint64_t source, uint64_t destination){
        EXPECT_EQ(source + destination,0);
    });
    //in range
    graph.rangeQuery(1,7, [](uint64_t timestamp,uint64_t source, uint64_t destination){
        EXPECT_GT(source + destination,0);
    });
}

//bigfile is reddit-dataset
TEST(FindTests, bigFile){
    //test size
    AdjList graph;
    graph.addFromFile("../../google_test/UnitTestInputs/largeFile");
    //expecting 5 timestamps
    EXPECT_EQ(graph.getSize(),5);
    //finds edge, edge exists
    EXPECT_TRUE(graph.findEdge(3626139277963222129, 4263427141073911189, 2014));
    //finds edge, edge does not exist
    EXPECT_FALSE(graph.findEdge(1,1,2014));
    //find edge, out of range
    EXPECT_FALSE(graph.findEdge(3626139277963222129, 4263427141073911189, 2017,2020));
    //find edge in range
    EXPECT_TRUE(graph.findEdge(3626139277963222129, 4263427141073911189, 2012,2020));

    //out of range
    graph.rangeQuery(2018,2020, [](uint64_t timestamp,uint64_t source, uint64_t destination){
        EXPECT_EQ(source + destination,0);
    });
    //in range
    graph.rangeQuery(2014,2018, [](uint64_t timestamp,uint64_t source, uint64_t destination){
        EXPECT_GT(source + destination,0);
    });
}

TEST(DeleteTest,bigFile){
//add and delete the whole file
    AdjList graph;
    graph.addFromFile("../../google_test/UnitTestInputs/largeFile");
    graph.addFromFile("../../google_test/UnitTestInputs/largeFileDel");
    EXPECT_EQ(graph.getSize(),0);
}

TEST(MemoryConsumption, bigFile){
    //test memory consumption of 1 Edge, 1 uint_64 takes 8 bytes; 1 edge has 24 bytes
    AdjList graph;
    graph.addFromFile("../../google_test/UnitTestInputs/largeFile");
    //count how many rows
    size_t result=0;
    for (int64_t i = 2013; i < 2018; ++i) {
        result+=graph.getEdgeCount(i);
    }
    EXPECT_EQ(graph.memoryConsumption(), result*24 );
}

