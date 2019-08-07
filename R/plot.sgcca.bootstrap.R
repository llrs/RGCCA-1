plot.bootstrap.sgcca = function(object, nb_bloc = length(object$A)){
  if (class(object) != "sgcca.bootstrap"){
    return("object is not sgcca.bootstrap !")
  }
  Stab = list()
  g = list()
  for (i in 1:nb_bloc){
    Stab[[i]] = data.frame(Variables = names(object$count[[i]])[1:object$top],Count = object$count[[i]][1:object$top]/object$nb_boot)
    g[[i]] = ggplot(Stab[[i]], aes(x = reorder(Variables, Count), y = Count))+
      coord_flip() +
      geom_bar(stat = "identity")+
      ggtitle(paste("Top" ,top, "bloc",i))+
      xlab("Variables")+
      ylab("Occurences")+
      theme(axis.title=element_text(size=15))+
      theme(axis.title.y=element_blank())+
      theme(plot.title=element_text(size=20))+
      theme(axis.text = element_text(size=15))
  }
  g
}